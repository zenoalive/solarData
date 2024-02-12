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
# Populate lists with planet names and their coordinates
planet_names = ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
longitude_values = [2.608642258043293, 0.5150643857021435, -2.843342599199367, -2.3633430084721843, 1.913159983202125, -0.9299284156186811, -3.080523873160146]
latitude_values = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

# Plot the positions of the planets
plt.scatter(longitude_values, latitude_values, label='Planets', color='blue', marker='o')

# Annotate each point with the corresponding planet name
for i, name in enumerate(planet_names):
    plt.annotate(name, (longitude_values[i], latitude_values[i]))

# Add a legend to the plot
plt.legend()

# Show the plot
# Define colors for each planet
colors = ['gray', 'orange', 'red', 'brown', 'gold', 'lightblue', 'blue']

# Plot the positions of the planets with colors
plt.figure(figsize=(10, 6))
for i, name in enumerate(planet_names):
    plt.scatter(longitude_values[i], latitude_values[i], label=name, color=colors[i], marker='o')

# Add labels and title
plt.xlabel('Geocentric Ecliptic Longitude (degrees)')
plt.ylabel('Geocentric Ecliptic Latitude (degrees)')
plt.title('Positions of Planets Relative to Earth')

# Add a legend to the plot
plt.legend()

# Show the plot
plt.grid(True)
# Plot the positions of the planets with colors and varying marker sizes based on distance
plt.figure(figsize=(10, 6))
for i, name in enumerate(planet_names):
    plt.scatter(longitude_values[i], latitude_values[i], label=name, color=colors[i], s=100 * distances_values[i], alpha=0.7)

# Add labels and title
plt.xlabel('Geocentric Ecliptic Longitude (degrees)')
plt.ylabel('Geocentric Ecliptic Latitude (degrees)')
plt.title('Positions of Planets Relative to Earth')

# Add a legend to the plot
plt.legend()

plt.show()
# plt.show()
