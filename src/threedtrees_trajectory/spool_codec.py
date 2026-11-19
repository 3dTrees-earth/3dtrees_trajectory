"""Lossless blocks for large tile spools with repeated byte patterns."""
from __future__ import annotations

import ctypes
import ctypes.util
import struct
import zlib
from functools import lru_cache

import numpy as np

HEADER = struct.Struct("<IIB")
ZSTANDARD, DEFLATE = 1, 2
MAX_BLOCK_BYTES = 256 * 1024**2


@lru_cache(maxsize=1)
def _zstandard():
    soname = ctypes.util.find_library("zstd")
    if soname is None:
        return None
    try:
        lib = ctypes.CDLL(soname)
    except OSError:
        return None
    size, pointer = ctypes.c_size_t, ctypes.c_void_p
    lib.ZSTD_compressBound.argtypes, lib.ZSTD_compressBound.restype = [size], size
    lib.ZSTD_compress.argtypes, lib.ZSTD_compress.restype = [pointer, size, pointer, size, ctypes.c_int], size
    lib.ZSTD_decompress.argtypes, lib.ZSTD_decompress.restype = [pointer, size, pointer, size], size
    lib.ZSTD_isError.argtypes, lib.ZSTD_isError.restype = [size], ctypes.c_uint
    return lib


def append_block(path, records):
    """Write one independently decodable byte-shuffled block."""
    raw = records.view(np.uint8).reshape(len(records), records.dtype.itemsize).T.copy()
    lib = _zstandard()
    if lib is None:
        data = zlib.compress(raw, level=1)
        method = DEFLATE
    else:
        encoded = np.empty(lib.ZSTD_compressBound(raw.nbytes), np.uint8)
        size = lib.ZSTD_compress(encoded.ctypes.data, encoded.nbytes, raw.ctypes.data, raw.nbytes, -1)
        if lib.ZSTD_isError(size):
            raise ValueError("Could not compress an occlusion spool block.")
        data = encoded[:size]
        method = ZSTANDARD
    with path.open("ab") as stream:
        stream.write(HEADER.pack(len(records), len(data), method))
        stream.write(data)
    return HEADER.size + len(data)


def _header(stream, width, end):
    raw = stream.read(HEADER.size)
    if not raw:
        return None
    if len(raw) != HEADER.size:
        raise ValueError("Incomplete occlusion spool block header.")
    count, size, method = HEADER.unpack(raw)
    if (not 0 < count * width <= MAX_BLOCK_BYTES or not 0 < size <= end - stream.tell()
            or size > MAX_BLOCK_BYTES + 4 * 1024**2):
        raise ValueError("Invalid occlusion spool block size.")
    return count, size, method


def block_count(path, dtype):
    count, width, end = 0, dtype.itemsize, path.stat().st_size
    with path.open("rb") as stream:
        while entry := _header(stream, width, end):
            rows, size, _ = entry
            count += rows
            stream.seek(size, 1)
    return count


def read_blocks(path, dtype):
    width, end = dtype.itemsize, path.stat().st_size
    with path.open("rb") as stream:
        while entry := _header(stream, width, end):
            count, size, method = entry
            encoded = stream.read(size)
            if method == ZSTANDARD:
                lib = _zstandard()
                if lib is None:
                    raise ValueError("Resuming this occlusion spool requires libzstd.")
                raw = np.empty(count * width, np.uint8)
                unpacked = lib.ZSTD_decompress(raw.ctypes.data, raw.nbytes, encoded, len(encoded))
                if lib.ZSTD_isError(unpacked) or unpacked != raw.nbytes:
                    raise ValueError("Could not decompress an occlusion spool block.")
            elif method == DEFLATE:
                try:
                    decoder = zlib.decompressobj()
                    raw = np.frombuffer(decoder.decompress(encoded, count * width + 1), np.uint8)
                except zlib.error as error:
                    raise ValueError("Could not decompress an occlusion spool block.") from error
                if not decoder.eof or len(raw) != count * width:
                    raise ValueError("Unexpected occlusion spool block length.")
            else:
                raise ValueError("Unknown occlusion spool block compression.")
            yield raw.reshape(width, count).T.copy().view(dtype).reshape(-1)
