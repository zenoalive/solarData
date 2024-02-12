
import ephem
import math
import time
import matplotlib.pyplot as plt
# Constants
G = 6.67430e-11  # Gravitational constant (m^3/kg/s^2)
M_sun = 1.989e30  # Mass of the Sun (kg)
# Define argument of perihelion values for each planet (in degrees)
argument_of_perihelion_values = {
    'Mercury': 77.45645,
    'Venus': 131.53298,
    'Mars': 336.04084,
    'Jupiter': 14.75385,
    'Saturn': 92.43194,
    'Uranus': 170.96424,
    'Neptune': 44.97135
}
colors = ['red', 'blue', 'green', 'orange', 'purple', 'yellow', 'cyan']


# Get heliocentric coordinates of earth
def calculate_earth_heliocentric_coordinates(observer):
    # Orbital elements of the Earth (taken from astronomical references)
    semi_major_axis = 1.000001018
    eccentricity = 0.016708634
    inclination = math.radians(0.00005)
    longitude_of_ascending_node = math.radians(-11.26064)
    argument_of_perihelion = math.radians(102.94719)
    mean_anomaly_at_epoch = math.radians(100.46435)

    # Julian date of epoch (J2000.0)
    epoch = 2451545.0

    # Current Julian date
    current_date = observer.date

    # Calculate time difference from epoch in days
    time_difference_days = current_date - epoch

    # Mean anomaly (M) at current date
    mean_anomaly = mean_anomaly_at_epoch + time_difference_days * (2 * math.pi / 365.256363)

    # Eccentric anomaly (E)
    eccentric_anomaly = solve_keplers_equation(mean_anomaly, eccentricity)

    # True anomaly (ν)
    true_anomaly = 2 * math.atan(math.sqrt((1 + eccentricity) / (1 - eccentricity)) * math.tan(eccentric_anomaly / 2))

    # Heliocentric distance (r)
    heliocentric_distance = semi_major_axis * (1 - eccentricity ** 2) / (1 + eccentricity * math.cos(true_anomaly))

    # Heliocentric ecliptic longitude (λ)
    heliocentric_longitude = longitude_of_ascending_node + argument_of_perihelion + true_anomaly

    return heliocentric_longitude, 0, heliocentric_distance  # Assuming Earth's heliocentric latitude is 0

# def calculate_earth_heliocentric_coordinates(observer):
#     # Create an Earth object
#     earth = ephem.EarthSatellite()

#     # Compute the position of the Earth
#     earth.compute(observer)

#     # Extract heliocentric ecliptic coordinates
#     heliocentric_longitude = earth.ecliptic_longitude
#     heliocentric_latitude = earth.ecliptic_latitude
#     heliocentric_distance = earth.sun_distance

#     return heliocentric_longitude, heliocentric_latitude, heliocentric_distance
# Function to get the argument of perihelion for a given planet
def get_argument_of_perihelion(planet_name):
    angles_in_radians = math.radians(argument_of_perihelion_values.get(planet_name, None))
    return angles_in_radians

# Example usage
# planet_name = 'Mercury'
# argument_of_perihelion = get_argument_of_perihelion(planet_name)
# print(argument_of_perihelion)  # Output: 77.45645


def get_mean_motion(planet):
    # Calculate orbital period (T) using Kepler's third law
    planet.compute()
    a = planet.sun_distance  # Semi-major axis (AU)
    T = 2 * math.pi * math.sqrt(a**3 / (G * M_sun))  # Orbital period (seconds)
    
    # Calculate mean motion (n)
    n = 2 * math.pi / T  # Mean motion (radians per second)
    
    return n

def get_mean_anomaly(planet, observer):
    planet.compute(observer)
    n = get_mean_motion(planet)  # Mean motion of the planet (angular speed)
    T0 = 2458426.053  # Approximate time of perihelion passage
    t = ephem.now()  # Current Julian date
    
    mean_anomaly = n * (t - T0)
    return mean_anomaly

def calculate_mean_anomalies(observer):
    planets = ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']
    mean_anomalies = {}
    
    for planet_name in planets:
        planet = getattr(ephem, planet_name)()
        mean_anomaly = get_mean_anomaly(planet, observer)
        mean_anomalies[planet_name] = mean_anomaly
    
    return mean_anomalies

def convert_to_geocentric(heliocentric_longitude, heliocentric_latitude, heliocentric_distance,
                          earth_heliocentric_longitude, earth_heliocentric_latitude, earth_heliocentric_distance):
    # Convert heliocentric coordinates to rectangular coordinates
    x_helio = heliocentric_distance * math.cos(heliocentric_latitude) * math.cos(heliocentric_longitude)
    y_helio = heliocentric_distance * math.cos(heliocentric_latitude) * math.sin(heliocentric_longitude)
    z_helio = heliocentric_distance * math.sin(heliocentric_latitude)

    # Convert Earth's heliocentric coordinates to rectangular coordinates
    x_earth = earth_heliocentric_distance * math.cos(earth_heliocentric_latitude) * math.cos(earth_heliocentric_longitude)
    y_earth = earth_heliocentric_distance * math.cos(earth_heliocentric_latitude) * math.sin(earth_heliocentric_longitude)
    z_earth = earth_heliocentric_distance * math.sin(earth_heliocentric_latitude)

    # Calculate geocentric coordinates
    x_geo = x_helio - x_earth
    y_geo = y_helio - y_earth
    z_geo = z_helio - z_earth

    # Convert rectangular coordinates back to spherical coordinates
    geocentric_distance = math.sqrt(x_geo ** 2 + y_geo ** 2 + z_geo ** 2)
    geocentric_longitude = math.atan2(y_geo, x_geo)
    geocentric_latitude = math.asin(z_geo / geocentric_distance)

    return geocentric_longitude, geocentric_latitude, geocentric_distance
# Eccentric anomaly
# def solve_keplers_equation(mean_anomaly, eccentricity):
#     # Initial guess for eccentric anomaly
#     E = mean_anomaly + eccentricity * math.sin(mean_anomaly)
    
#     # Iterate to solve Kepler's equation
#     max_iterations = 1000
#     tolerance = 1e-10
#     for _ in range(max_iterations):
#         next_E = E - (E - eccentricity * math.sin(E) - mean_anomaly) / (1 - eccentricity * math.cos(E))
#         if abs(next_E - E) < tolerance:
#             break
#         E = next_E
    
#     return E
def solve_keplers_equation(mean_anomaly, eccentricity):
    # Initial guess for eccentric anomaly
    E = mean_anomaly
    
    # Iterate to solve Kepler's equation
    max_iterations = 1000
    tolerance = 1e-10
    for _ in range(max_iterations):
        next_E = E - (E - eccentricity * math.sin(E) - mean_anomaly) / (1 - eccentricity * math.cos(E))
        if abs(next_E - E) < tolerance:
            break
        E = next_E
    
    return E



def calculate_eccentric_anomaly(planet, observer):
    planet.compute(observer)
    mean_anomaly = get_mean_anomaly(planet, observer)
    eccentricity = 1 - planet.sun_distance / planet.earth_distance #eccentricity
    
    # Solve Kepler's equation for eccentric anomaly
    eccentric_anomaly = solve_keplers_equation(mean_anomaly, eccentricity)
    
    return eccentric_anomaly

# True Anomaly

def calculate_true_anomaly(eccentricity, eccentric_anomaly):
    # Calculate true anomaly (ν) using Kepler's equation
    tan_half_nu = math.sqrt((1 + eccentricity) / (1 - eccentricity)) * math.tan(eccentric_anomaly / 2)
    nu = 2 * math.atan(tan_half_nu)
    
    return nu

# Calculate the distance
def calculate_distance(semi_major_axis, eccentricity, true_anomaly):
    # Calculate the distance (r) using the formula
    r = semi_major_axis * (1 - eccentricity**2) / (1 + eccentricity * math.cos(true_anomaly))
    
    return r

def calculate_ecliptic_longitude(true_anomaly, argument_of_perihelion):
    
    # Calculate the heliocentric ecliptic longitude (λ)
    longitude = true_anomaly + argument_of_perihelion
    return longitude

def calculate_ecliptic_latitude():
    # Since the orbits lie in the ecliptic plane, the latitude (β) is always zero
    return 0.0



if __name__ == "__main__":
    observer = ephem.Observer()
    observer.lat = '51.5074'  # Latitude of the observer (e.g., London)
    observer.lon = '0.1278'   # Longitude of the observer (e.g., London)
    observer.elevation = 10    # Elevation of the observer (e.g., in meters)
    observer.date = '2024/02/12 12:00:00'  # Date and time of observation

# Calculate Earth's heliocentric coordinates
    earth_heliocentric_longitude, earth_heliocentric_latitude, earth_heliocentric_distance = calculate_earth_heliocentric_coordinates(observer)

# Print the results
    print("Heliocentric Ecliptic Longitude of Earth:", earth_heliocentric_longitude)
    print("Heliocentric Ecliptic Latitude of Earth:", earth_heliocentric_latitude)
    print("Distance from Sun to Earth:", earth_heliocentric_distance)

    # Calculate mean anomalies
    mean_anomalies = calculate_mean_anomalies(observer)
    planets = ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']

    # Calculate eccentric anomalies
    eccentric_anomalies = {}
    for planet_name in planets:
        planet = getattr(ephem, planet_name)()
        eccentric_anomalies[planet_name] = calculate_eccentric_anomaly(planet, observer)
    

    
    #calculate True anomalies
    true_anomalies = {}
    heliocentric_distances = {}
    heliocentric_longitudes = {}
    heliocentric_latitudes = {}
    for planet_name in planets :
        planet = getattr(ephem, planet_name)()
        eccentric_anomaly = calculate_eccentric_anomaly(planet, observer)  # Calculate eccentric anomaly first
        eccentricity = 1 - planet.sun_distance / planet.earth_distance  # Calculate eccentricity
        true_anomaly = calculate_true_anomaly(eccentricity, eccentric_anomaly)  # Calculate true anomaly
        true_anomalies[planet_name] = true_anomaly
        w = get_argument_of_perihelion(planet_name)
        a = planet.sun_distance
        heliocentric_distances[planet_name] = calculate_distance(a,eccentricity, true_anomaly )
        heliocentric_longitudes[planet_name] = calculate_ecliptic_longitude(true_anomaly, w)
        heliocentric_latitudes[planet_name] = calculate_ecliptic_latitude()    
    
    # Plot the positions of the planets with colors and varying marker sizes based on distance
    plt.figure(figsize=(10, 6))
    for i, name in enumerate(planets):
        plt.scatter(heliocentric_longitudes[name], heliocentric_latitudes[name], label=name, color=colors[i], s=100 * heliocentric_distances[name], alpha=0.7)

    # Add labels and title
    plt.xlabel('Geocentric Ecliptic Longitude (degrees)')
    plt.ylabel('Geocentric Ecliptic Latitude (degrees)')
    plt.title('Positions of Planets Relative to Earth')

# Add a legend to the plot
    plt.legend()
    plt.show()

   


    # Print mean anomalies and eccentric anomalies for comparison
    # for planet_name in planets:
    #     print(f"{planet_name}:")
        # print(f"Mean Anomaly: {mean_anomalies[planet_name]}")
        # print(f"Eccentric Anomaly: {eccentric_anomalies[planet_name]}")
        # print(f"True Anomaly: {true_anomalies[planet_name]}")
        # print(f"Distance from sun: {distances[planet_name]}")
        # print(f"Heliocentric ecliptic longitude is : {heliocentric_longitudes[planet_name]}")
        # print()


## Geocentric part
    geocentric_coordinates = {}

    for planet_name in planets:
        heliocentric_longitude = heliocentric_longitudes[planet_name]
        heliocentric_latitude = heliocentric_latitudes[planet_name]
        heliocentric_distance = heliocentric_distances[planet_name]

    # Call the function to convert to geocentric coordinates
        geocentric_longitude, geocentric_latitude, geocentric_distance = convert_to_geocentric(
        heliocentric_longitude, heliocentric_latitude, heliocentric_distance,
        earth_heliocentric_longitude, earth_heliocentric_latitude, earth_heliocentric_distance
        )

    # Store the geocentric coordinates
        geocentric_coordinates[planet_name] = (geocentric_longitude, geocentric_latitude, geocentric_distance)

    # Print the geocentric coordinates
    for planet_name, coordinates in geocentric_coordinates.items():
        print(f"{planet_name}:")
        print(f"Geocentric Ecliptic Longitude: {coordinates[0]}")
        print(f"Geocentric Ecliptic Latitude: {coordinates[1]}")
        print(f"Distance from Earth: {coordinates[2]}")
        print()
print('It worked')
