# DeepEntropyX QTVI pipeline

**Confidential. Proprietary to DEMAZUMDER LLC. Do not distribute.**

A Qt desktop tool for marking noise and cardiac features in polysomnography and
bedside-monitor recordings. It reads EDF or `.dat` recordings, converts them to a
packed `.bin` format, and presents them for manual annotation across ECG, PPG,
arterial-pressure, accelerometer, respiratory, EEG and hypnogram channels. A
second stage then marks per-beat templates.

Supported datasets: **MESA**, **BITTIUM**, **CHAOS**.

---

## Running the prebuilt executable

1. Download `exe_files.zip` and unzip it.
2. Run `noise_marking_gui.exe` **from inside that folder** — it needs the Qt
   `.dll` files and `config.csv` that sit alongside it.
3. When prompted, type `1` for MESA, `2` for BITTIUM, or `3` for CHAOS.

---

## Building from source

### Requirements

|                  |                                                                                  |
| ---------------- | -------------------------------------------------------------------------------- |
| **C++ standard** | **C++20** (`/std:c++20` on MSVC, `-std=c++20` on GCC/Clang)                      |
| **Compiler**     | MSVC 2022 (v143 toolset). GCC 11+ or Clang 14+ should also work but are untested |
| **Qt**           | **6.7 or newer**, modules `Core` `Gui` `Widgets` `Charts` `OpenGLWidgets`        |
| **Build system** | Visual Studio + Qt VS Tools, or CMake 3.19+                                      |

Two notes on those versions, since both are easy to get wrong:

- **Qt Charts is a separate installer component.** It is not part of the base Qt
  install. If a fresh clone fails with a missing `QtCharts/QChartView`, re-run
  the Qt Maintenance Tool and tick *Qt Charts* under your Qt 6.x version.
- **`OpenGLWidgets` is its own module in Qt 6.** `QOpenGLWidget` moved out of
  `Widgets`, so linking only `Qt6::Widgets` produces unresolved symbols.

The Qt 6.7 floor comes from the `.ui` files, which are saved with fully-qualified
enum names (`Qt::Orientation::Horizontal`). Older versions of `uic` cannot parse
them.

### Bundled third-party code

These live in the source tree; nothing needs to be fetched separately.

- `file_format_parsing/edflib.h` — EDF/BDF reading (C, compile as C or with
  `extern "C"`)
- `file_format_parsing/pugixml.hpp` — XML parsing for sleep-stage sidecars

### Build steps — Visual Studio

1. Install the **Qt VS Tools** extension.
2. *Extensions → Qt VS Tools → Qt Versions* → add your Qt 6.7+ `msvc2022_64`
   installation.
3. Open the solution, then *Project → Properties*:
   - *C/C++ → Language → C++ Language Standard* → `/std:c++20`
   - *Qt Project Settings → Qt Modules* → `core gui widgets charts openglwidgets`
   - *Debugging → Working Directory* → `$(OutDir)`, so that `config.csv` is
     found at runtime (see the note below)
4. Build. Qt VS Tools runs `uic` on the `.ui` files to generate
   `ui_noise_marking_gui.h` and `ui_TemplateViewerWindow.h`; do not commit or
   hand-edit those.

### Build steps — CMake

```cmake
cmake_minimum_required(VERSION 3.19)
project(noise_marking_gui LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOUIC ON)      # generates ui_*.h from the .ui files
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 6.7 REQUIRED COMPONENTS Core Gui Widgets Charts OpenGLWidgets)

qt_add_executable(noise_marking_gui
    main.cpp gui_handler.cpp signal_renderer.cpp
    user_annotation_handler.cpp user_control_handler.cpp user_marking_handler.cpp
    annotation_eraser.cpp grid_overlay.cpp bin_chunk_loader.cpp
    config_loader.cpp file_to_bin.cpp anneal_handler.cpp
    BinPlotWidget.cpp FocusPanelWidget.cpp
    noise_marking_gui.ui
)

qt_add_executable(template_marking
    template_marking.cpp TemplateViewerWindow.cpp feature_marks.cpp
    config_loader.cpp
    TemplateViewerWindow.ui
)

foreach(tgt noise_marking_gui template_marking)
    target_link_libraries(${tgt} PRIVATE
        Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Charts Qt6::OpenGLWidgets)
    target_include_directories(${tgt} PRIVATE ${CMAKE_SOURCE_DIR})
endforeach()
```

Configure and build:

```bash
cmake -B build -DCMAKE_PREFIX_PATH="C:/Qt/6.7.2/msvc2022_64"
cmake --build build --config Release
```

Adjust the source lists to match your tree — the headers expect the
subdirectories `file_format_parsing/`, `logging/`, `peak_finding/`, `theme/` and
`template_marking_gui/` to be reachable from the include path.

---

## Pointing the application at your data

All paths come from **`config.csv`**, one row per dataset, selected by the
`data_type` column matching your `1`/`2`/`3` choice at startup.

`config.csv` is opened as a **relative path**, so it must be in the working
directory the program is launched from — not necessarily next to the
executable. If you see `ERROR: cannot open config.csv`, that is the cause. In
Visual Studio the debugger defaults to `$(ProjectDir)`, so either set the working
directory to `$(OutDir)` or copy `config.csv` into the project root.

### The two columns you need

| Column               | Meaning                                              |
| -------------------- | ---------------------------------------------------- |
| `original_file_path` | Folder holding the source `.edf` / `.dat` recordings |
| `output_folder`      | Root folder for everything the tool writes           |

**Leave either blank and you will be prompted with a folder picker at startup.**
Filling them in with absolute paths skips the dialogs — the faster option if you
run repeatedly against the same dataset.

### Output layout

Subfolders are created automatically under `output_folder`:

```
output_folder/
├── noise_marking_output/     annotations from stage 1
├── r_peak_finding_output/    peak-finder results
├── annealed_output/
├── template_outputs/         *_templates.bin, consumed by template_marking
├── qtvi_marker_path/         template marker positions
├── quality_metric/
├── training_log/
└── snapshot_path/            saved plot images
```

### Other config columns

- `main_file_extention`, `sleep_file_extention` — which files to scan for
- `bin_size_minutes` — length of each processing chunk
- `*_raw_rate` / `*_upsampled_rate` — per-channel native and target sample rates.
  Every channel is resampled to its own target rate; a rate of `0` marks the
  channel absent.
- `notch_filter_hz` — powerline notch. `0` disables the toggle entirely.
- `waveform_highpass_hz`, `blanking_period`, `threshold` — peak-finder defaults
- `use_consensus_rpeak` — accepts `1`/`true`/`yes` or `0`/`false`/`no`
- `sleepstate_length` — sleep-stage epoch length in seconds (30s in MESA)
- `age`, `sex`, `height_cm`, `weight_kg`, `hr_rest`, `hr_max` — these will be added to later and are currently stored but for now, you can leave them blank

Channel *labels* are not in the CSV. They are hardcoded per dataset in
`apply_dataset_specific_channel_labels()` in `config_loader.cpp` — for example
MESA's ECG channel is `EKG` while BITTIUM's is `ECG_1`. Adding a new dataset
means editing that function as well as adding a CSV row.

---

## Annotating

### Marking noise

There are two ways to create an annotation.

**Drag and release** over the region you want to mark. This annotates one channel
at a time by default. To widen the scope, set the radio button to *Mark all Chan*
for every channel, or *Mark all ECG* for the ECG channels only. 

**Click 'Create Annotation'** to place separate start and stop markers, following
whatever the radio button is set to. Because the two clicks are independent, you
can set the start marker, scroll, and then set the end marker — useful for
regions longer than one screen.

Right-click a marker to remove it.

### Navigation

`←` or `A` scrolls back, `→` or `D` scrolls forward, by the step in the skip
interval box.

### The grid

The *grid* button overlays thick lines every 0.2 s and thin lines every 0.04 s,
rescaling with the visible window length.

### Peak finder, blanking and threshold

The built-in peak identifier is **not** the full algorithm, so poor results here
do not mean the real peak finder will do badly.

To improve the visible GUI peak finding, Click the
*blanking and threshold* button, drag over a region, and set each value between
0 and 1 in the popup.

- **Blanking** is the interval after an R peak during which no further R peak may
  be detected. This is a raw ms value - it does not depend on previous segments of the dataset.
- **Threshold** is the height an R peak must reach. This is relative to the previous stretch of data.

Marking noise helps here too: regions marked as noise are excluded when blanking and threshold are estimated.

### Saving

Click *save* when noise and feature marking are done. Processing starts
immediately — some of it must finish before you continue, the rest runs in
parallel while you work on the next step. Expect roughly a minute.

---

## Template marking

After processing, a new window opens for template marking, with nine markers.
If they are too crowded, un-check the PPG or ECG marker boxes to hide a set.

To flag a bad beat: **right-click once** for a bad R, **right-click twice** for a
bad PPG.

| #   | Marker  | Meaning                 |
| --- | ------- | ----------------------- |
| 1   | P beg   | Onset of P wave         |
| 2   | P peak  | Peak of P wave          |
| 3   | Q beg   | Onset of Q wave         |
| 4   | S end   | End of S wave (J Point) |
| 5   | T begin | Before the T wave       |
| 6   | T end   | After the T wave        |
| 7   | On      | PPG onset               |
| 8   | Pk      | PPG Systolic peak       |
| 9   | Dc      | Dicrotic notch          |
| 10  | En      | End                     |

Save when finished and the next file opens automatically.