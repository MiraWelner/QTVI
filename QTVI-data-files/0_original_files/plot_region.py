import matplotlib.pyplot as plt
import pandas as pd


def plot_detrended_data(file_path):
    # Calculate the starting index: 256 Hz * 490 minutes * 60 seconds
    start_index = 256 * 501 * 60
    end_index = 256 * 503 * 60

    try:
        # Load the CSV file.
        # Since it's likely a signal file, we assume it contains a single column without a header.
        # If the file has a header, change header=None to header=0.
        data = pd.read_csv(file_path, header=None)

        # Extract values from the start index to the end
        plot_data = data.iloc[start_index:end_index]

        # Check if we have data to plot
        if plot_data.empty:
            print(f"No data found starting from index {start_index}.")
            return

        # Plot the data
        plt.figure(figsize=(24, 6))
        plt.plot(plot_data.index, plot_data.values, label="Detrended Signal")
        plt.title(f"Detrended Data starting from index {start_index}")
        plt.xlabel("Index (256Hz)")
        plt.ylabel("Value")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    csv_file = "3010201_20120320_annealed_detrended.csv"
    plot_detrended_data(csv_file)
    plot_detrended_data("3010201_20120320_threshold_0.40_detrended.csv")
