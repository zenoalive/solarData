import ephem
import math
import datetime
import matplotlib.pyplot as plt

# Define the list of planet names
planet_names = ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
planet_colors = {'Mercury': 'gray', 'Venus': 'orange', 'Mars': 'red', 'Jupiter': 'brown', 'Saturn': 'gold', 'Uranus': 'cyan', 'Neptune': 'blue'}

# Initialize a dictionary to store latitudes and longitudes for each planet
planet_positions = {name: {'latitudes': [], 'longitudes': []} for name in planet_names}

# Define observer
observer = ephem.Observer()

# Define the start and end dates for the year 2023
start_date = datetime.date(2023, 1, 1)
end_date = datetime.date(2023, 12, 31)
dates = [start_date + datetime.timedelta(days=i) for i in range((end_date - start_date).days + 1)]

# Define the time interval as one day
delta = datetime.timedelta(days=1)

# Iterate over each day in the year 2023
current_date = start_date

while current_date <= end_date:
    # Set the observer date
    observer.date = current_date
    
    # Iterate over each planet
    for planet_name in planet_names:
        # Calculate the position of the current planet
        planet = getattr(ephem, planet_name)()
        planet.compute(observer)
        
        # Store the latitude and longitude in the dictionary
        planet_positions[planet_name]['latitudes'].append(planet.hlat)
        planet_positions[planet_name]['longitudes'].append(planet.hlong)
    
    # Move to the next day
    current_date += delta
    
    # Plot the orbits for each planet
    plt.figure(figsize=(8, 6))
    for planet_name in planet_names:
        plt.scatter(planet_positions[planet_name]['latitudes'], planet_positions[planet_name]['longitudes'], marker='.', color=planet_colors[planet_name], label=planet_name)
    
    # Add labels and legend
    plt.xlabel('Heliocentric Ecliptic Latitude (degrees)')
    plt.ylabel('Heliocentric Ecliptic Longitude (degrees)')
    plt.title('Orbital Paths of Planets in 2023')
    plt.legend()
    
    # Show the plot
    plt.show()
    
    # Close the current figure
    plt.close()

