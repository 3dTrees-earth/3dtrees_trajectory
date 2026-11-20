from setuptools import Extension, setup

setup(ext_modules=[
    Extension(
        "threedtrees_trajectory._ray_walk",
        ["src/threedtrees_trajectory/ray_walk.cpp"],
        language="c++",
        extra_compile_args=["-std=c++17", "-O3", "-fno-fast-math", "-ffp-contract=off"],
    ),
])
