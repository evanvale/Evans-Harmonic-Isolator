# Evan's Harmonic Isolator

A CLAP audio plugin that isolates harmonic content with tunable frequency selection and musical intervals.

## Features

- **Root Note Selection**: A3 to C6 with fine tuning (±100 cents)
- **Musical Tuning System**: Semitone-based frequency selection
- **Harmonic Content**: Isolate 5th, 9th, and 13th harmonics
- **Octave Controls**: Independent low and high octave adjustment
- **Audio Processing**: 2x oversampling, saturation, stereo spread
- **Performance**: SIMD optimized with parameter smoothing

## Parameters

| Parameter | Range | Description |
|-----------|--------|-------------|
| **Root Note (440Hz)** | A3-C6 | Base frequency for harmonic isolation |
| **Fine Tune** | ±100 cents | Fine pitch adjustment in cents |
| **Bandwidth** | 1-45% | Filter bandwidth (inverse Q factor) |
| **Root Gain** | 0-100% | Amplitude of the root frequency |
| **Harmonics** | 0-100% | Amount of musical harmonics (5th, 9th, 13th) |
| **Low Octaves** | 0-100% | Sub-octave content below root |
| **High Octaves** | 0-100% | Super-octave content above root |
| **Spread** | L100%-R100% | Stereo spread amount |
| **Saturation** | 0-10x | Harmonic saturation drive |
| **LPF Cutoff** | 40Hz-20kHz | Low-pass filter frequency |
| **Wet Boost** | -10 to +20dB | Output level compensation |
| **Dry/Wet** | 0-100% | Mix between processed and dry signal |

## Installation

Download the `.clap` file for your platform from [Releases](../../releases) and copy to your CLAP plugin folder.

## Platform Support

Builds automatically for:
- **Linux**: x86_64
- **Windows**: x86_64  
- **macOS**: Intel (x86_64) and Apple Silicon (ARM64)

## License

Public domain. Do whatever you want with this.

## Credits

- Development collaboration between Evan and Claude AI
- [CLAP](https://github.com/free-audio/clap) - CLever Audio Plugin API
- GitHub Actions for automated builds
