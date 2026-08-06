# %% imports
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]


# ===================================================
# Helper functions
# ===================================================
def circular_mean(angles):
    # 1.convert to radians
    rads = np.radians(angles)

    # 2. Duplicate angle to account for circular nature
    rads_x2 = rads * 2

    # 3. Vector components
    mean_sin = np.mean(np.sin(rads_x2))
    mean_cos = np.mean(np.cos(rads_x2))

    # 4. Get mean angle back in rads
    mean_rad = np.arctan2(mean_sin, mean_cos) / 2

    # 5. Get mean angle back in degrees
    return np.degrees(mean_rad)


def frames_to_minutes(num_frames):
    time_minutes = np.zeros(num_frames)
    for i in range(num_frames):
        if i <= 300:
            time_sec = i * 2
        else:
            time_sec = 600 + (i - 300) * 10
        time_minutes[i] = time_sec / 60.0
    return time_minutes


def time_based_moving_average(time_sec, data, window_sec=15.0):
    """
    Applies a rolling window average over the data over a specific time window.
    NOTE: This is currently ONLY used for smoothing out the Trajectory plots
    (Angle SD vs Max Count) to avoid chaotic visual squiggles.
    """
    smoothed = np.zeros(len(data))
    for i, t in enumerate(time_sec):
        valid_idx = (time_sec >= t - window_sec / 2) & (time_sec <= t + window_sec / 2)
        window_data = data[valid_idx]
        if np.all(np.isnan(window_data)):
            smoothed[i] = np.nan
        else:
            smoothed[i] = np.nanmean(window_data)
    return smoothed


def circular_std(angles):
    """
    Calculates the Circular Standard Deviation of a list of angles.
    Because angles wrap around (0 degrees = 180 degrees for axial data like actin fibers),
    we use vector math to find the resultant length R, and map it to a standard deviation.
    A low SD means fibers are parallel (ordered). A high SD means fibers are crossed or swirled (chaotic).
    """
    # 1. Convert to radians
    rads = np.radians(angles)

    # 2. Duplicate angle to account for axial nature (0=180)
    rads_x2 = rads * 2

    # 3. Vector components to find resultant vector length R
    mean_sin = np.mean(np.sin(rads_x2))
    mean_cos = np.mean(np.cos(rads_x2))
    R = np.sqrt(mean_sin**2 + mean_cos**2)

    # 4. Calculate circular standard deviation
    # bounded to avoid division by zero or log of 0
    if R <= 0.0:
        return np.inf
    std_rad = (
        np.sqrt(-np.log(R)) / 2
    )  # von misses distribution approximation for axial data

    # 5. Get std back in degrees
    return np.degrees(std_rad)


def extract_metadata(folder_name):
    # splits the folder string
    parts = folder_name.split("_")
    if len(parts) >= 2:
        experiment_lipid = parts[0]
        frog = parts[1]

        # join the rest of the strings in case there are more underscores
        condition = "_".join(parts[2:]) if len(parts) > 2 else ""
        return experiment_lipid, frog, condition
    else:
        return None, None, None


# =============================================
# Folder processing
# =============================================
def process_folder(folder_path, divisions_per_axis=2):
    """
    Reads all CSV files from a specific folder (representing an experiment movie).
    Groups data by frame (Slice) and subdivides the space into 'divisions_per_axis'^2 regions.
    Extracts the dominant orientation (angle) and strength (max count) of actin fibers
    for each sub-region, filtering out noise using a Coherency threshold.
    Returns arrays of Angle Standard Deviation (to detect Swirls) and Average Max Count.
    """
    # divisions_per_axis = 2 creates a 2x2 grid (4 sub-regions) INSIDE each CSV.
    # With 4 CSV files, this results in 16 total regions across the whole image.
    # Change to 4 for a 4x4 grid per CSV (64 total), or 8 for an 8x8 grid (256 total).

    csv_list = list(folder_path.glob("*.csv"))

    if len(csv_list) == 0:
        return None, None

    slice_data = {}

    for csv_idx, file in enumerate(csv_list):
        df = pd.read_csv(file)

        # --- EXCEPTION FOR MOVIES THAT NEED CUTS ---
        if "D_FROG153-2_rGDB-GFP750nm_Utr594100nm_rep1" in folder_path.name:
            df = df[df["Slice"] >= 242]
        # --- DYNAMIC SUB-QUADRANT CREATION ---
        # Divides the specific CSV (which is already 1 quadrant) into smaller bins
        df["X_bin"] = pd.cut(df["X"], bins=divisions_per_axis, labels=False)
        df["Y_bin"] = pd.cut(df["Y"], bins=divisions_per_axis, labels=False)

        # Calculate a local ID for the sub-region within this specific CSV (e.g., 1 to 4)
        local_pq = (df["Y_bin"] * divisions_per_axis + df["X_bin"]) + 1

        # Calculate a Global PQ ID to prevent overlap between the 4 CSV files
        # Example for divisions=2: CSV 0 gets IDs 1-4, CSV 1 gets 5-8, CSV 2 gets 9-12...
        regions_per_csv = divisions_per_axis**2
        df["Global_PQ"] = (csv_idx * regions_per_csv) + local_pq

        # Group data by temporal Frame (Slice) and the unique Global Pseudo-Quadrant
        for key, group in df.groupby(["Slice", "Global_PQ"]):
            slice_val = key[0] if isinstance(key, tuple) else key
            slice_num = int(float(str(slice_val)))

            if slice_num not in slice_data:
                slice_data[slice_num] = {
                    "angles": [],
                    "counts": [],
                    "total_quads": 0,
                }

            slice_data[slice_num]["total_quads"] += 1
            mean_coherency = group["Coherency"].mean()

            # --- ANTI-NOISE FILTER ---
            if mean_coherency >= 0.055:
                angles_mod = np.mod(group["Orientation"], 180)

                y, _ = np.histogram(angles_mod, bins=np.arange(-2.5, 182.6, 5))
                y[0] += y[-1]
                y = y[:-1]
                x_centers = np.arange(0, 180, 5)

                max_count = np.max(y)
                ymax_idx = np.argmax(y)
                most_angle = x_centers[ymax_idx]

                slice_data[slice_num]["angles"].append(most_angle)
                slice_data[slice_num]["counts"].append(max_count)
            else:
                slice_data[slice_num]["counts"].append(0)

    if not slice_data:
        return None, None

    num_frames = max(slice_data.keys()) + 1
    sd_angle = np.zeros(num_frames)
    average_count = np.zeros(num_frames)

    for i in range(num_frames):
        if i in slice_data:
            data = slice_data[i]

            if data["total_quads"] > 0:
                average_count[i] = sum(data["counts"]) / data["total_quads"]
            else:
                average_count[i] = 0

            # Calculate Circular SD using all valid sub-regions across all 4 CSVs
            # Require at least 4 valid sub-regions to compute a meaningful SD
            if len(data["angles"]) >= 4:
                sd_angle[i] = circular_std(data["angles"])
            else:
                sd_angle[i] = np.nan
        else:
            sd_angle[i] = np.nan
            average_count[i] = np.nan

    return sd_angle, average_count


# ======================================================================
# Graphing and main function
# ======================================================================
def four_cornered_analysis():
    """
    Main execution pipeline.
    1. Scans the directory for experiment folders.
    2. Processes each folder to extract Angle SD and Max Count over time.
    3. Groups the data by Lipid Condition.
    4. Generates 2x2 grid plots for timeline comparisons.
    5. Generates a Trajectory plot (SD vs Max Count) to visualize organization dynamics.
    """
    main_folder = Path(".")
    experiment_folders = set(archivo.parent for archivo in main_folder.rglob("*.csv"))
    # Group data by variable
    data_per_variable = {}
    # Group per defined variables
    for folder in experiment_folders:
        lipid, frog, condition = extract_metadata(folder.name)

        # skip if not matching the expected pattern
        if frog is None or lipid is None:
            print(f"skipping folder {folder} due to unexpected naming pattern.")
            continue
        sd_angle, average_count = process_folder(folder)

        if average_count is None or sd_angle is None:
            print(f"Skipping folder {folder} due to no CSV files.")
            continue
        if average_count is not None and sd_angle is not None:
            if lipid not in data_per_variable:
                data_per_variable[lipid] = []
            data_per_variable[lipid].append((frog, condition, sd_angle, average_count))
    # 2. Plotting Setup
    image_groups = [
        ("A_SM", ["A", "SM"]),
        ("F_D", ["F", "D"]),
    ]

    cmap = plt.get_cmap("tab10")
    # Map each unique frog base name to a consistent color
    unique_frogs = set()
    for lipid in data_per_variable:
        for frog, _, _, _ in data_per_variable[lipid]:
            unique_frogs.add(frog)
    unique_frogs = sorted(list(unique_frogs))
    frog_colors = {frog: cmap(i % 10) for i, frog in enumerate(unique_frogs)}

    # ======================================================================
    # 3. Draw the Timeline Graphs (Time vs Max Count & Time vs Angle SD)
    # ======================================================================
    for group_name, target_lipids in image_groups:
        fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharex=True)

        for col_idx, current_lipid in enumerate(target_lipids):
            ax_count = axes[0, col_idx]
            ax_peak = axes[1, col_idx]

            ax_count.set_title(f"Lipid Condition: {current_lipid}", fontsize=26)
            ax_count.set_ylabel("Average Max Count", fontsize=22)
            ax_count.set_ylim(0, 700)
            ax_count.tick_params(axis="both", labelsize=20)
            ax_count.grid(True, linestyle="--", alpha=0.5)

            ax_peak.set_xlabel("Time (Minutes)", fontsize=22)
            ax_peak.set_ylabel("Angle SD (Degrees) - Swirl/Chaos", fontsize=22)
            ax_peak.set_ylim(0, 80)
            ax_peak.tick_params(axis="both", labelsize=20)
            ax_peak.grid(True, linestyle="--", alpha=0.5)

            if current_lipid in data_per_variable:
                for line_idx, (frog, condition, sd_angle, avg_count) in enumerate(
                    data_per_variable[current_lipid]
                ):
                    label_name = frog
                    linestyle = "-"  # Default to solid line

                    if "REP" in condition:
                        rep_parts = [p for p in condition.split("_") if "REP" in p]
                        if rep_parts:
                            label_name += f" ({rep_parts[-1]})"
                            if rep_parts[-1] != "REP1":
                                linestyle = "--"

                    time_minutes = frames_to_minutes(len(avg_count))
                    color = frog_colors[frog]

                    # Plot Max Count
                    ax_count.plot(
                        time_minutes,
                        avg_count,
                        linestyle=linestyle,
                        color=color,
                        marker="o",
                        markersize=5,
                        alpha=0.8,
                        label=label_name,
                    )

                    # Plot Angle SD (Swirled vs Ordered)
                    # Plot the line for everything (matplotlib auto-clips at y=80 limit)
                    ax_peak.plot(
                        time_minutes,
                        sd_angle,
                        linestyle=linestyle,
                        color=color,
                        alpha=0.8,
                        label=label_name,
                    )

                    # Plot points only for valid values <= 80 to avoid clipping artifacts
                    valid_idx = np.where(~np.isnan(sd_angle) & (sd_angle <= 80))[0]
                    if len(valid_idx) > 0:
                        ax_peak.plot(
                            time_minutes[valid_idx],
                            sd_angle[valid_idx],
                            linestyle="None",
                            color=color,
                            marker="o",
                            markersize=5,
                            alpha=0.8,
                        )

        # Master Legend for the Figure
        lines_labels = [ax.get_legend_handles_labels() for ax in axes.flatten()]
        lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
        unique_legend = dict(zip(labels, lines))
        if unique_legend:
            fig.legend(
                unique_legend.values(),
                unique_legend.keys(),
                loc="lower center",
                fontsize=16,
                title_fontsize=16,
                bbox_to_anchor=(0.5, 0.0),
                ncol=min(len(unique_legend), 4),
                columnspacing=1.0,
                title="Experiments",
            )


        plt.tight_layout(rect=(0, 0.12, 1.0, 0.96))
        plt.savefig(f"Lipid_Experiments_{group_name}.png", dpi=300)

    # ======================================================================
    # 4. Draw the Angle SD vs Max Count correlation graphs (Trajectory)
    # [WORK IN PROGRESS] - This visualization logic is still experimental
    # ======================================================================
    fig_corr, axes_corr = plt.subplots(2, 2, figsize=(16, 12))
    all_lipids = ["A", "SM", "F", "D"]

    for idx, current_lipid in enumerate(all_lipids):
        ax = axes_corr.flatten()[idx]
        ax.set_title(f"Lipid Condition: {current_lipid}", fontsize=26)
        ax.set_xlabel("Average Max Count", fontsize=22)
        ax.set_ylabel("Angle SD (Degrees) - Swirl/Chaos", fontsize=22)
        ax.set_ylim(0, 80)
        ax.set_xlim(0, 700)
        ax.tick_params(axis="both", labelsize=20)
        ax.grid(True, linestyle="--", alpha=0.5)

        if current_lipid in data_per_variable:
            for frog, condition, sd_angle, avg_count in data_per_variable[
                current_lipid
            ]:
                label_name = frog
                if "REP" in condition:
                    rep_parts = [p for p in condition.split("_") if "REP" in p]
                    if rep_parts:
                        label_name += f" ({rep_parts[-1]})"

                color = frog_colors[frog]
                time_minutes = frames_to_minutes(len(avg_count))

                # Smooth the data for the trajectory plot to avoid chaotic squiggles
                smooth_sd = time_based_moving_average(
                    time_minutes, sd_angle, window_sec=30.0
                )
                smooth_count = time_based_moving_average(
                    time_minutes, avg_count, window_sec=30.0
                )

                # Draw the trajectory line
                ax.plot(
                    smooth_count,
                    smooth_sd,
                    linestyle="-",
                    color=color,
                    alpha=0.6,
                    label=label_name,
                )

                # Mark the Start (S) and End (E) of the trajectory
                valid_indices = np.where(
                    ~np.isnan(smooth_sd) & ~np.isnan(smooth_count)
                )[0]
                if len(valid_indices) > 0:
                    start_idx = valid_indices[0]
                    end_idx = valid_indices[-1]

                    # Draw Start marker (Circle)
                    ax.plot(
                        smooth_count[start_idx],
                        smooth_sd[start_idx],
                        marker="o",
                        color=color,
                        markersize=8,
                        alpha=0.9,
                        markeredgecolor="black",
                    )
                    ax.text(
                        smooth_count[start_idx],
                        smooth_sd[start_idx],
                        " S",
                        fontsize=16,
                        color="black",
                        weight="bold",
                    )

                    # Draw End marker (Square)
                    ax.plot(
                        smooth_count[end_idx],
                        smooth_sd[end_idx],
                        marker="s",
                        color=color,
                        markersize=8,
                        alpha=0.9,
                        markeredgecolor="black",
                    )
                    ax.text(
                        smooth_count[end_idx],
                        smooth_sd[end_idx],
                        " E",
                        fontsize=16,
                        color="black",
                        weight="bold",
                    )

    # Add master legend for trajectory figure
    lines_labels = [ax.get_legend_handles_labels() for ax in axes_corr.flatten()]
    lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
    unique_legend = dict(zip(labels, lines))
    if unique_legend:
        fig_corr.legend(
            unique_legend.values(),
            unique_legend.keys(),
            loc="lower center",
            fontsize=16,
            title_fontsize=16,
            bbox_to_anchor=(0.5, 0.0),
            ncol=min(len(unique_legend), 4),
            columnspacing=1.0,
            title="Experiments",
        )


    plt.tight_layout(rect=(0, 0.12, 1.0, 0.96))
    plt.savefig("Lipid_Experiments_Trajectory.png", dpi=300)


if __name__ == "__main__":
    four_cornered_analysis()
