import matplotlib.pyplot as plt

# Define the size of the plot
plt.figure(figsize=(10, 6))

# Set the title and labels
plt.title("Geocentric Positions of Planets (Date: YYYY-MM-DD)")
plt.xlabel("Geocentric Ecliptic Longitude (degrees)")
plt.ylabel("Geocentric Ecliptic Latitude (degrees)")

# Enable gridlines for better visualization
plt.grid(True)

# Initialize lists to store planet names and their coordinates
planet_names = []
longitude_values = []
latitude_values = []
