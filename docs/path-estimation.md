# Path estimation

[Back to README](../README.md)

| Input evidence | Method |
| --- | --- |
| Enough multi-return pulses with GPS time | `lines` fits position and velocity in each window using pulse endpoints. Pairing uses time, point source and scanner channel within each file. |
| Too many unpaired LAS/LAZ returns | The reader tries a 1 µs time tolerance for complete, unique return sequences on the same source and channel. |
| Usable scan angles and GPS time | `scan_angle` fits the angle off nadir, assuming a level line scanner. If pairing failed, the affected file goes to the single-return traversal. |
| Missing path windows in a file with usable single returns | `lines_scan_angle` keeps reliable line fits and uses scan angles to fill missing time windows. The combined estimate is marked `rough`. |
| A delivered path with timestamps, in any EPSG CRS or local like the cloud | Occlusion uses the delivered positions. An estimated trajectory is optional. |
| Too little geometric information for either estimator | The refusal message asks for the delivered trajectory. |

Line fits that miss their pulse lines by more than ten times the median, with a threshold of at least 0.25 m, are
dropped. The tool also drops positions over 1 m from the path described by nearby fits. The line quality labels are
`good` through a 0.1 m median miss and `rough` through 0.5 m. A line estimate above 0.5 m falls back to scan angles if available. Scan-angle fits are labelled `rough`. On the archived 3Dtrees comparison, most YellowScan multi-return paths were accurate to a few
centimetres. Scan-angle estimates had about 1 m horizontal error on RIEGL and DJI flights and about 8 m on
YellowScan.

The tool checks for static scans using the largest block by file, point source and time. It looks for a single
position that explains the pulse lines within 5 cm. For clouds without enough pulse lines, it checks whether
near-ground points concentrate around one location at least 15 times more than an even spread would.
`--no-platform-check` is a CLI diagnostic for known moving acquisitions. Static scan position estimation is outside
the path tool's scope.

## Pulse pairing

A pulse is the set of echoes of one laser shot: the same GPS time, point source and scanner channel, with return
numbers 1 to NumberOfReturns. The reader pairs consecutive records in scanner order. Echoes of one shot that are
stored apart are joined after reading by GPS time, point source and channel within each file ("stray" echoes, up to
40 million per file). When too many echoes stay unpaired, the reader retries with a 1 µs time tolerance.

Some exports keep only the **first and the declared last return** of a shot and drop intermediate returns, and also
store them apart (dataset 2587, LiDAR360: groups of return 1 and 3 with NumberOfReturns 3, or a lone return 3). When
strict pairing fails for such a file, the separated echoes are joined again allowing missing intermediate returns
(return 1 and the declared last return present, no repeated return number). Those first-to-last lines feed the path
estimate; the file is still traced as single-return for occlusion. Files that pair normally are unchanged. On 2587
this gives the `lines` method with 1.7 mm median distance to the delivered path.

A file whose pairing fails and that has neither such lines nor a usable scan angle cannot be estimated (exit 4)
unless a delivered path is given.

## Quality

| Label | Meaning |
| --- | --- |
| `good` | `lines`, median miss of the pulse lines to the fitted path ≤ 0.1 m |
| `rough` | `lines` with median miss 0.1–0.5 m, or any scan-angle position |
| `unreliable` | median miss > 0.5 m; the tool retries with scan angles when the file has them |

The **miss** is the perpendicular distance from the fitted sensor position to the line through the first and last
return of a pulse. Accuracy against delivered paths: [validation.md](validation.md).

## Static check

The tool checks for static scans using the largest block by file, point source and time. It looks for a single
position that explains the pulse lines within 5 cm. For clouds without enough pulse lines, it checks whether
near-ground points concentrate around one location at least 15 times more than an even spread would.
`--no-platform-check` is a CLI diagnostic for known moving acquisitions. Static scan position estimation is outside
the path tool's scope.
