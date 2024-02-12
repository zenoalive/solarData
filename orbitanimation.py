import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import ephem
import datetime

# Function to update the plot at each frame

planet_names = ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
planet_colors = {'Mercury': 'gray', 'Venus': 'orange', 'Mars': 'red', 'Jupiter': 'brown', 'Saturn': 'gold', 'Uranus': 'cyan', 'Neptune': 'blue'}
# Define observer
observer = ephem.Observer()

# Define the start and end dates for the year 2023
start_date = datetime.date(2023, 1, 1)
end_date = datetime.date(2023, 12, 31)
dates = [start_date + datetime.timedelta(days=i) for i in range((end_date - start_date).days + 1)]
# Calculate the total number of days between start_date and end_date
total_days = (end_date - start_date).days

# Assuming a frame rate of 30 frames per second (fps)
frame_rate = 30

# Calculate the number of frames
num_frames = total_days * frame_rate

# Define the time interval as one day
delta = datetime.timedelta(days=1)

# Iterate over each day in the year 2023
current_date = start_date
# Convert the datetime.timedelta to an ephem.Date
delta_days = delta.days  # Extract the number of days from the timedelta
delta_ephem = ephem.Date(delta_days)  # Convert the number of days to an ephem.Date object

# Update the observer date

def update(frame):
    # Clear the previous plot
    plt.clf()
    
    # Set the observer date
    observer.date += delta_ephem
    
    # Plot the orbits for each planet
    for planet_name in planet_names:
        planet = getattr(ephem, planet_name)()
        planet.compute(observer)
        plt.scatter(planet.hlat, planet.hlong, marker='.', color=planet_colors[planet_name], label=planet_name, s=300)
    
    # Add labels and legend
    plt.xlabel('Heliocentric Ecliptic Latitude (degrees)')
    plt.ylabel('Heliocentric Ecliptic Longitude (degrees)')
    plt.title('Orbital Paths of Planets in 2023')
    plt.legend()

# Initialize the observer and set the start date
observer = ephem.Observer()
observer.date = start_date

# Initialize the figure
fig = plt.figure(figsize=(8, 6))

# Create an animation
ani = FuncAnimation(fig, update, frames=num_frames, interval=100)
ani.save('solar_system_animation.gif', writer='imagemagick')

# Show the animation
plt.show()
# 