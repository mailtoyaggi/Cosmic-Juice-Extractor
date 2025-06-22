import ephem
from math import degrees, fabs
from datetime import datetime, timezone, timedelta
from timezonefinder import TimezoneFinder
import requests
import swisseph as swe
import pytz


# Function to calculate aspects
def calculate_aspect(angle1, angle2, is_sun_or_moon):
    angle_diff = fabs(angle1 - angle2)
    if angle_diff > 180:
        angle_diff = 360 - angle_diff

    aspects = {
        "Conjunction": 0,
        "Opposition": 180,
        "Trine": 120,
        "Square": 90,
        "Sextile": 60,
        "Quincunx": 150,
        "Semi-Square": 45,
        "Sesquiquadrate": 135,
        "Semi-Sextile": 30,
    }

    orbs = {
        "Conjunction": 8 if is_sun_or_moon else 6,
        "Opposition": 8 if is_sun_or_moon else 6,
        "Trine": 8 if is_sun_or_moon else 6,
        "Square": 8 if is_sun_or_moon else 6,
        "Sextile": 4,
        "Quincunx": 3,
        "Semi-Square": 3,
        "Sesquiquadrate": 3,
        "Semi-Sextile": 2,
    }

    closest_aspect = min(aspects,
                         key=lambda asp: fabs(aspects[asp] - angle_diff))
    orb = fabs(aspects[closest_aspect] - angle_diff)

    if orb <= orbs[closest_aspect]:
        return closest_aspect, orb
    else:
        return "No major aspect", orb

def adjust_degree(degree, adjustment):
    # Subtract adjustment and ensure the degree wraps around within 0 to 360
    adjusted_degree = (degree - adjustment) % 360
    return adjusted_degree

def get_zodiac(degree):
    # Adjust degree before categorizing
    adjusted_degree = adjust_degree(degree, 0)
    signs = [
        ('Aries', 0, 29),
        ('Taurus', 30, 59),
        ('Gemini', 60, 89),
        ('Cancer', 90, 119),
        ('Leo', 120, 149),
        ('Virgo', 150, 179),
        ('Libra', 180, 209),
        ('Scorpio', 210, 239),
        ('Sagittarius', 240, 269),
        ('Capricorn', 270, 299),
        ('Aquarius', 300, 329),
        ('Pisces', 330, 359)
    ]

    for sign, start, end in signs:
        if start <= adjusted_degree <= end:
            degree_in_sign = adjusted_degree - start
            return f"{sign} {degree_in_sign:.2f}°"
    return "Unknown"  # In case the degree is outside 0–359°

def get_nakshatra(degree):
    # Adjust degree before categorizing
    adjusted_degree = adjust_degree(degree, 0)
    nakshatras = [
        ('Ashwini', 0, 13.33),
        ('Bharani', 13.33, 26.67),
        ('Krittika', 26.67, 40.00),
        ('Rohini', 40.00, 53.33),
        ('Mrigashira', 53.33, 66.67),
        ('Ardra', 66.67, 80.00),
        ('Punarvasu', 80.00, 93.33),
        ('Pushya', 93.33, 106.67),
        ('Ashlesha', 106.67, 120.00),
        ('Magha', 120.00, 133.33),
        ('Purva Phalguni', 133.33, 146.67),
        ('Uttara Phalguni', 146.67, 160.00),
        ('Hasta', 160.00, 173.33),
        ('Chitra', 173.33, 186.67),
        ('Swati', 186.67, 200.00),
        ('Vishakha', 200.00, 213.33),
        ('Anuradha', 213.33, 226.67),
        ('Jyeshtha', 226.67, 240.00),
        ('Mula', 240.00, 253.33),
        ('Purva Ashadha', 253.33, 266.67),
        ('Uttara Ashadha', 266.67, 280.00),
        ('Shravana', 280.00, 293.33),
        ('Dhanishta', 293.33, 306.67),
        ('Shatabhisha', 306.67, 320.00),
        ('Purva Bhadrapada', 320.00, 333.33),
        ('Uttara Bhadrapada', 333.33, 346.67),
        ('Revati', 346.67, 360.00)
    ]

    for nakshatra, start, end in nakshatras:
        if start <= adjusted_degree < end:
            degree_in_nakshatra = adjusted_degree - start
            return f"{nakshatra} {degree_in_nakshatra:.2f}°"
    return "Unknown"  # In case the degree is outside 0–360°

def calculate_house_cusps(ascendant):
    """
    Calculate house cusps starting from the Ascendant degree.

    Args:
        ascendant (float): The Ascendant degree in sidereal mode.

    Returns:
        dict: A dictionary of house numbers (1–12) and their cusp degrees.
    """
    house_cusps = {}
    for house_number in range(1, 13):
        cusp_degree = (ascendant + (house_number - 1) * 30) % 360  # Add 30° for each house and wrap around 360°
        house_cusps[house_number] = cusp_degree
    return house_cusps

def find_house_for_planet(planet_degree, house_cusps):
    """
    Determine which house a planet falls into and how many degrees it is into that house.

    Args:
        planet_degree (float): The degree position of the planet.
        house_cusps (dict): A dictionary of house numbers and their cusp degrees.

    Returns:
        tuple: The house number and degrees into that house.
    """
    sorted_houses = sorted(house_cusps.items(), key=lambda x: x[1])  # Sort houses by cusp degree
    for i, (house, cusp) in enumerate(sorted_houses):
        next_cusp = sorted_houses[(i + 1) % 12][1]  # Wrap around to the first cusp after the 12th house
        if cusp <= planet_degree < next_cusp or (cusp > next_cusp and (planet_degree >= cusp or planet_degree < next_cusp)):
            degrees_into_house = (planet_degree - cusp) % 360  # Calculate how far the planet is into the house
            return house, degrees_into_house

def get_ascendant_midheaven(date_str, time_str, latitude, longitude):
    # Extract UTC date and time
    year, month, day = map(int, date_str.split('/'))
    hour, minute, second = map(int, time_str.split(':'))
    hour = hour + minute / 60 + second / 3600  # Convert to decimal hours

    # Convert to Julian Day Number
    julian_day = swe.julday(year, month, day, hour)

    # Set Lahiri Ayanamsa for sidereal calculations
    swe.set_sid_mode(swe.SIDM_LAHIRI)

    # Calculate houses & angles
    houses, angles = swe.houses(julian_day, latitude, longitude, b'P')  # Placidus system

    # Extract Ascendant and Midheaven
    ascendant = angles[0]  # Ascendant (ASC)
    midheaven = angles[1]  # Midheaven (MC)

    # Get the Lahiri ayanamsa for sidereal adjustment
    ayanamsa = swe.get_ayanamsa(julian_day)

    # Adjust Ascendant and Midheaven to sidereal mode
    ascendant_sidereal = (ascendant - ayanamsa) % 360
    midheaven_sidereal = (midheaven - ayanamsa) % 360

    return ascendant_sidereal, midheaven_sidereal

# Function to compare aspects within natal planetary positions
def calculate_and_sort_natal_aspects(positions):
    """Calculate aspects between natal planetary positions and organize them per planet."""
    natal_aspects = {planet: [] for planet in positions.keys()}  # Initialize dictionary

    planets = list(positions.keys())
    for i in range(len(planets)):
        for j in range(i + 1, len(planets)):
            planet1, degree1 = planets[i], positions[planets[i]]
            planet2, degree2 = planets[j], positions[planets[j]]
            aspect, orb = calculate_aspect(degree1, degree2, False)

            if aspect != "No major aspect":
                natal_aspects[planet1].append(f"{aspect} {planet2} (Orb: {orb:.2f}°)")
                natal_aspects[planet2].append(f"{aspect} {planet1} (Orb: {orb:.2f}°)")

    return natal_aspects

def calculate_part_of_fortune(ascendant, sun, moon, is_day_chart):
    """
    Calculate the Part of Fortune.

    Args:
        ascendant (float): The Ascendant degree in sidereal mode.
        sun (float): The Sun's degree in sidereal mode.
        moon (float): The Moon's degree in sidereal mode.
        is_day_chart (bool): True if it's a day chart, False if it's a night chart.

    Returns:
        float: The sidereal degree of the Part of Fortune.
    """
    if is_day_chart:
        part_of_fortune = ascendant + moon - sun
    else:
        part_of_fortune = ascendant - moon + sun

    # Ensure the degree is within 0–360°
    return part_of_fortune % 360


def get_vimshottari_dashas(moon_degree, birth_datetime):
    """Return Vimshottari Maha Dashas and Antar Dashas from birth for 120 years."""

    nakshatras = [
        ('Ashwini', 0), ('Bharani', 13.33), ('Krittika', 26.67), ('Rohini', 40.00),
        ('Mrigashira', 53.33), ('Ardra', 66.67), ('Punarvasu', 80.00), ('Pushya', 93.33),
        ('Ashlesha', 106.67), ('Magha', 120.00), ('Purva Phalguni', 133.33),
        ('Uttara Phalguni', 146.67), ('Hasta', 160.00), ('Chitra', 173.33),
        ('Swati', 186.67), ('Vishakha', 200.00), ('Anuradha', 213.33), ('Jyeshtha', 226.67),
        ('Mula', 240.00), ('Purva Ashadha', 253.33), ('Uttara Ashadha', 266.67),
        ('Shravana', 280.00), ('Dhanishta', 293.33), ('Shatabhisha', 306.67),
        ('Purva Bhadrapada', 320.00), ('Uttara Bhadrapada', 333.33), ('Revati', 346.67)
    ]

    dasha_sequence = [
        ('Ketu', 7), ('Venus', 20), ('Sun', 6), ('Moon', 10), ('Mars', 7),
        ('Rahu', 18), ('Jupiter', 16), ('Saturn', 19), ('Mercury', 17)
    ]
    dasha_years = dict(dasha_sequence)
    nakshatra_rulers = [seq[0] for seq in dasha_sequence] * 3  # To map 27 nakshatras

    # Find nakshatra index and start degree
    for idx, (_, start_deg) in enumerate(nakshatras):
        end_deg = nakshatras[(idx + 1) % len(nakshatras)][1]
        if start_deg <= moon_degree < end_deg or (idx == len(nakshatras) - 1 and moon_degree >= start_deg):
            nakshatra_index = idx
            nakshatra_start = start_deg
            break

    ruler = nakshatra_rulers[nakshatra_index]
    ruler_index = next(i for i, (p, _) in enumerate(dasha_sequence) if p == ruler)

    degrees_into_nakshatra = moon_degree - nakshatra_start
    portion_elapsed = degrees_into_nakshatra / 13.33
    balance_portion = 1 - portion_elapsed
    total_years = dasha_years[ruler]
    balance_years = total_years * balance_portion

    dashas = []
    current_date = birth_datetime

    for i in range(len(dasha_sequence)):
        maha_planet, maha_years = dasha_sequence[(ruler_index + i) % len(dasha_sequence)]

        if i == 0:
            dasha_span = timedelta(days=balance_years * 365.25)
            full_span_years = maha_years
            maha_elapsed_portion = portion_elapsed
        else:
            dasha_span = timedelta(days=maha_years * 365.25)
            full_span_years = maha_years
            maha_elapsed_portion = 0

        dasha_start = current_date
        dasha_end = dasha_start + dasha_span

        # Rotate antar sequence to start from maha lord
        maha_index = next(j for j, (p, _) in enumerate(dasha_sequence) if p == maha_planet)
        rotated_sequence = dasha_sequence[maha_index:] + dasha_sequence[:maha_index]

        antars = []
        antar_start = dasha_start
        days_elapsed = maha_elapsed_portion * maha_years * 365.25 if i == 0 else 0
        skipped = 0

        for antar_planet, antar_years in rotated_sequence:
            antar_fraction = antar_years / 120.0
            antar_days = antar_fraction * full_span_years * 365.25

            if i == 0 and skipped + antar_days < days_elapsed:
                skipped += antar_days
                continue

            if i == 0 and skipped < days_elapsed:
                leftover = antar_days - (days_elapsed - skipped)
                antar_span = timedelta(days=leftover)
                antar_start = dasha_start  # First antar starts at actual birth time
            else:
                antar_span = timedelta(days=antar_days)

            antar_end = antar_start + antar_span
            if antar_start < dasha_end:
                trimmed_end = min(antar_end, dasha_end)
                antars.append((antar_planet, antar_start, trimmed_end))
                antar_start = trimmed_end
            else:
                break

        dashas.append((maha_planet, dasha_start, dasha_end, antars))
        current_date = dasha_end

    return dashas


def estimate_orbital_period_days(planet_id):
    # Sidereal orbital periods in days (approximate)
    sidereal_days = {
        swe.MOON: 27.32166,
        swe.SUN: 365.25636,
        swe.MERCURY: 87.969,
        swe.VENUS: 224.701,
        swe.MARS: 686.98,
        swe.JUPITER: 4332.59,
        swe.SATURN: 10759.22,
        swe.URANUS: 30685.4,
        swe.NEPTUNE: 60189.0,
        swe.PLUTO: 90560.0,
        swe.MEAN_NODE: 6798.38,  # Rahu/Ketu
        swe.TRUE_NODE: 6798.38
    }
    return sidereal_days.get(planet_id, 365.25)


def calculate_return_age(planet_id, birth_longitude, birth_julian, end_julian, orb=0.01, cooldown_deg=10.0, step_days=1.0):
    """
    Counts the number of times a planet has returned to its natal longitude between two Julian dates.
    - Uses sidereal positions from Swiss Ephemeris.
    - Accounts for retrogrades using cooldown angle.
    """
    swe.set_sid_mode(swe.SIDM_LAHIRI)

    count = 0
    last_return_deg = None
    current_jd = birth_julian + step_days

    while current_jd <= end_julian:
        pos, _ = swe.calc_ut(current_jd, planet_id, swe.FLG_SIDEREAL)
        current_long = pos[0]
        delta = abs((current_long - birth_longitude + 360) % 360)

        if delta <= orb:
            if last_return_deg is None or abs((current_long - last_return_deg + 360) % 360) > cooldown_deg:
                count += 1
                last_return_deg = current_long

        current_jd += step_days

    # Compute fractional return for current progress
    total_days = end_julian - birth_julian
    estimated_period = estimate_orbital_period_days(planet_id)
    return round(total_days / estimated_period, 2)



def determine_lordships(ascendant_degree):
    """Determine planetary lordships based on the Ascendant sign."""
    ascendant_sign = get_zodiac(ascendant_degree).split()[0]  # Extract sign name

    # Define rulerships for each Ascendant sign
    rulerships = {
        "Aries": {1: "Mars", 2: "Venus", 3: "Mercury", 4: "Moon", 5: "Sun", 6: "Mercury",
                  7: "Venus", 8: "Mars", 9: "Jupiter", 10: "Saturn", 11: "Saturn", 12: "Jupiter"},
        "Taurus": {1: "Venus", 2: "Mercury", 3: "Moon", 4: "Sun", 5: "Mercury", 6: "Venus",
                   7: "Mars", 8: "Jupiter", 9: "Saturn", 10: "Saturn", 11: "Jupiter", 12: "Mars"},
        "Gemini": {1: "Mercury", 2: "Moon", 3: "Sun", 4: "Mercury", 5: "Venus", 6: "Mars",
                   7: "Jupiter", 8: "Saturn", 9: "Saturn", 10: "Jupiter", 11: "Mars", 12: "Venus"},
        "Cancer": {1: "Moon", 2: "Sun", 3: "Mercury", 4: "Moon", 5: "Mars", 6: "Jupiter",
                   7: "Saturn", 8: "Saturn", 9: "Jupiter", 10: "Mars", 11: "Venus", 12: "Mercury"},
        "Leo": {1: "Sun", 2: "Mercury", 3: "Venus", 4: "Mars", 5: "Jupiter", 6: "Saturn",
                7: "Saturn", 8: "Jupiter", 9: "Mars", 10: "Venus", 11: "Mercury", 12: "Moon"},
        "Virgo": {1: "Mercury", 2: "Venus", 3: "Mars", 4: "Jupiter", 5: "Saturn", 6: "Saturn",
                  7: "Jupiter", 8: "Mars", 9: "Venus", 10: "Mercury", 11: "Moon", 12: "Sun"},
        "Libra": {1: "Venus", 2: "Mars", 3: "Jupiter", 4: "Saturn", 5: "Saturn", 6: "Jupiter",
                  7: "Mars", 8: "Venus", 9: "Mercury", 10: "Moon", 11: "Sun", 12: "Mercury"},
        "Scorpio": {1: "Mars", 2: "Jupiter", 3: "Saturn", 4: "Saturn", 5: "Jupiter", 6: "Mars",
                    7: "Venus", 8: "Mercury", 9: "Moon", 10: "Sun", 11: "Mercury", 12: "Venus"},
        "Sagittarius": {1: "Jupiter", 2: "Saturn", 3: "Saturn", 4: "Jupiter", 5: "Mars", 6: "Venus",
                        7: "Mercury", 8: "Moon", 9: "Sun", 10: "Mercury", 11: "Venus", 12: "Mars"},
        "Capricorn": {1: "Saturn", 2: "Saturn", 3: "Jupiter", 4: "Mars", 5: "Venus", 6: "Mercury",
                      7: "Moon", 8: "Sun", 9: "Mercury", 10: "Venus", 11: "Mars", 12: "Jupiter"},
        "Aquarius": {1: "Saturn", 2: "Jupiter", 3: "Mars", 4: "Venus", 5: "Mercury", 6: "Moon",
                     7: "Sun", 8: "Mercury", 9: "Venus", 10: "Mars", 11: "Jupiter", 12: "Saturn"},
        "Pisces": {1: "Jupiter", 2: "Mars", 3: "Venus", 4: "Mercury", 5: "Moon", 6: "Sun",
                   7: "Mercury", 8: "Venus", 9: "Mars", 10: "Jupiter", 11: "Saturn", 12: "Saturn"},
    }

    return rulerships.get(ascendant_sign, {})

# Function to get planetary positions
def get_planetary_positions(date_str, time_str, latitude, longitude, include_asc_midheaven=False):
    # Convert date/time to Julian Day
    year, month, day = map(int, date_str.split('/'))
    hour, minute, second = map(int, time_str.split(':'))
    hour_decimal = hour + minute / 60 + second / 3600
    julian_day = swe.julday(year, month, day, hour_decimal)

    # Set Lahiri Ayanamsa for sidereal calculations
    swe.set_sid_mode(swe.SIDM_LAHIRI)

    # Calculate sidereal planetary positions
    planets = {
        'Sun': swe.calc(julian_day, swe.SUN, swe.FLG_SIDEREAL),
        'Moon': swe.calc(julian_day, swe.MOON, swe.FLG_SIDEREAL),
        'Mercury': swe.calc(julian_day, swe.MERCURY, swe.FLG_SIDEREAL),
        'Venus': swe.calc(julian_day, swe.VENUS, swe.FLG_SIDEREAL),
        'Mars': swe.calc(julian_day, swe.MARS, swe.FLG_SIDEREAL),
        'Jupiter': swe.calc(julian_day, swe.JUPITER, swe.FLG_SIDEREAL),
        'Saturn': swe.calc(julian_day, swe.SATURN, swe.FLG_SIDEREAL),
        'Uranus': swe.calc(julian_day, swe.URANUS, swe.FLG_SIDEREAL),
        'Neptune': swe.calc(julian_day, swe.NEPTUNE, swe.FLG_SIDEREAL),
        'Pluto': swe.calc(julian_day, swe.PLUTO, swe.FLG_SIDEREAL),
        'Rahu': swe.calc(julian_day, swe.TRUE_NODE, swe.FLG_SIDEREAL),  # True Lunar Node
    }

    # Calculate Ketu (180° opposite Rahu)
    rahu_position = planets['Rahu'][0][0]  # Sidereal longitude of Rahu
    ketu_position = (rahu_position + 180) % 360  # Wrap within 0–360°
    planets['Ketu'] = ((ketu_position,),)  # Format similarly to other planetary data

    # Extract sidereal longitudes
    raw_degrees = {
        planet: planet_data[0][0]
        for planet, planet_data in planets.items()
    }

    # Calculate Ascendant and add to raw_degrees
    houses, angles = swe.houses_ex(julian_day, float(latitude), float(longitude), b'P', swe.FLG_SIDEREAL)
    ascendant_sidereal = angles[0]  # Ascendant (House 1 cusp)
    raw_degrees['Ascendant'] = ascendant_sidereal

    # Determine if it's a day chart (Sun above the horizon)
    is_day_chart = raw_degrees['Sun'] >= raw_degrees['Ascendant'] and raw_degrees['Sun'] < ((raw_degrees['Ascendant'] + 180) % 360)

    # Calculate the Part of Fortune
    part_of_fortune_degree = calculate_part_of_fortune(
        raw_degrees['Ascendant'], raw_degrees['Sun'], raw_degrees['Moon'], is_day_chart)
    raw_degrees['Part of Fortune'] = part_of_fortune_degree

    # Format degrees with zodiac and nakshatra information
    formatted_degrees = {
        planet: f"{raw_degrees[planet]:.2f}° ({get_zodiac(raw_degrees[planet])}, {get_nakshatra(raw_degrees[planet])})"
        for planet in raw_degrees
    }

    # Include Ascendant and Midheaven if requested
    if include_asc_midheaven:
        midheaven_sidereal = angles[1]  # Midheaven (House 10 cusp)
        raw_degrees['Midheaven'] = midheaven_sidereal
        formatted_degrees['Ascendant'] = f"{ascendant_sidereal:.2f}° ({get_zodiac(ascendant_sidereal)}, {get_nakshatra(ascendant_sidereal)})"
        formatted_degrees['Midheaven'] = f"{midheaven_sidereal:.2f}° ({get_zodiac(midheaven_sidereal)}, {get_nakshatra(midheaven_sidereal)})"

    # Determine planetary lordships based on Ascendant
    lordships = determine_lordships(ascendant_sidereal)

    # Append lordship information to formatted degrees
    for house, ruler in lordships.items():
        if ruler in formatted_degrees:
            formatted_degrees[ruler] += f" (Lord of {house}th house)"

    # Calculate aspects within the natal chart
    natal_aspects = calculate_and_sort_natal_aspects(raw_degrees)

    return raw_degrees, formatted_degrees, natal_aspects


# Function to convert local time to UTC
def convert_local_to_utc(date_str, time_str, latitude, longitude):
    tf = TimezoneFinder()
    timezone_name = tf.timezone_at(lat=float(latitude), lng=float(longitude))
    if not timezone_name:
        raise ValueError("Could not determine time zone for the given latitude and longitude.")

    from pytz import timezone as pytz_timezone
    local_timezone = pytz_timezone(timezone_name)

    # Handle dashes from HTML input and pad seconds if missing
    date_str = date_str.replace("-", "/")
    if len(time_str.strip().split(":")) == 2:
        time_str += ":00"

    local_time = datetime.strptime(date_str + ' ' + time_str, '%Y/%m/%d %H:%M:%S')
    local_time = local_timezone.localize(local_time)
    utc_time = local_time.astimezone(timezone.utc)

    return utc_time.strftime('%Y/%m/%d'), utc_time.strftime('%H:%M:%S')

def map_planets_to_houses(planets, house_cusps):
    """
    Determines which house each planet falls into and how many degrees into the house it is.
    Args:
        planets (dict): Dictionary of planet names and their degrees.
        house_cusps (dict): Dictionary of house numbers and their cusp degrees.
    Returns:
        dict: A dictionary mapping planets to the houses and degrees into the house.
    """
    planet_to_house = {}
    for planet, degree in planets.items():
        for house, cusp in house_cusps.items():
            next_cusp = house_cusps[(house % 12) + 1]  # Wrap around after 12th house
            if cusp <= degree < next_cusp or (cusp > next_cusp and (degree >= cusp or degree < next_cusp)):
                degrees_into_house = (degree - cusp) % 360  # Calculate how far the planet is into the house
                planet_to_house[planet] = (house, degrees_into_house)
                break
    return planet_to_house

def calculate_composite_chart(raw_positions1, raw_positions2):
    """
    Calculates the composite chart by finding the midpoints of planetary positions,
    Ascendant, and house cusps for two individuals.

    Args:
        raw_positions1 (dict): Dictionary of planetary positions for Person A.
        raw_positions2 (dict): Dictionary of planetary positions for Person B.

    Returns:
        tuple: A dictionary of composite planetary positions and composite house cusps.
    """
    composite_positions = {}

    # Calculate midpoints for planets and angles
    for planet in raw_positions1.keys():
        degree_a = raw_positions1[planet]
        degree_b = raw_positions2[planet]
        midpoint = (degree_a + degree_b) / 2
        if fabs(degree_a - degree_b) > 180:  # Handle 360° wrap-around
            midpoint = (degree_a + degree_b + 360) / 2 % 360
        composite_positions[planet] = midpoint

    # Calculate composite Ascendant as House 1 cusp
    composite_ascendant = composite_positions['Ascendant']

    # Generate composite house cusps starting from composite Ascendant
    composite_house_cusps = calculate_house_cusps(composite_ascendant)

    return composite_positions, composite_house_cusps

# Function to send a message via Telegram
def send_telegram_message(token, chat_id, message):
    url = f"https://api.telegram.org/bot{token}/sendMessage"
    payload = {'chat_id': chat_id, 'text': message}
    response = requests.post(url, data=payload)
    if response.status_code != 200:
        print(f"Failed to send message: {response.text}")
    else:
        print("Message sent successfully.")

def send_telegram_file(token, chat_id, file_content):
    url = f"https://api.telegram.org/bot{token}/sendDocument"
    files = {'document': ('astrology_report.txt', file_content.getvalue())}  # File name and content
    payload = {'chat_id': chat_id}
    response = requests.post(url, data=payload, files=files)
    if response.status_code != 200:
        print(f"Failed to send file: {response.text}")
    else:
        print("File sent successfully.")

# Function to compare planetary positions
def compare_planetary_positions(raw_positions1, raw_positions2):
    aspects = []
    for planet1, degree1 in raw_positions1.items():
        for planet2, degree2 in raw_positions2.items():
            aspect, orb = calculate_aspect(degree1, degree2, False)
            if aspect != "No major aspect":
                aspects.append((planet1, planet2, aspect, orb))
    # Sort aspects by orb (ascending)
    aspects.sort(key=lambda x: x[3])  # 4th element is the orb
    return aspects

import io

def main():
    # Default birth details in IST (India Standard Time)
    default_birth_date = '1996/06/02'
    default_birth_time = '18:26:00'
    default_latitude = '12.9716'
    default_longitude = '77.5946'

    # Ask user for input
    print("Do you want to use your stored birth details or enter new details?")
    print("1. Use my stored details")
    print("2. Enter new details")
    choice = input("Enter your choice (1 or 2): ")

    if choice == '1':
        birth_date = default_birth_date
        birth_time = default_birth_time
        latitude = default_latitude
        longitude = default_longitude
    elif choice == '2':
        birth_date = input("Enter birth date in YYYY/MM/DD format: ")
        birth_time = input("Enter birth time in HH:MM:SS format: ")
        latitude = input("Enter latitude (decimal format, e.g., 12.9716): ")
        longitude = input("Enter longitude (decimal format, e.g., 77.5946): ")
    else:
        print("Invalid choice! Exiting the program.")
        return

    print("Do you want to compare the chosen birth chart with today's date or another custom date?")
    print("1. Compare with today's date")
    print("2. Compare with a custom date")
    date_choice = input("Enter your choice (1 or 2): ")

    if date_choice == '1':
        today_ist = datetime.now().strftime('%Y/%m/%d')
        fixed_time_ist = '12:00:00'
        custom_latitude = latitude
        custom_longitude = longitude
    elif date_choice == '2':
        custom_date = input("Enter the custom date in YYYY/MM/DD format: ")
        try:
            # Validate the provided custom date
            custom_date_obj = datetime.strptime(custom_date, '%Y/%m/%d')
            today_ist = custom_date
            fixed_time_ist = input(
                "Enter the custom time in HH:MM:SS format (default is 12:00:00): "
            ) or '12:00:00'

            # Ask for latitude and longitude inputs for the custom date
            custom_latitude = input("Enter latitude for custom date (decimal format, e.g., 12.9716): ")
            custom_longitude = input("Enter longitude for custom date (decimal format, e.g., 77.5946): ")
        except ValueError:
            print("Invalid date format! Please enter the date in YYYY/MM/DD format.")
            return
    else:
        print("Invalid choice! Exiting the program.")
        return

    print("Thank you! Processing your input...")

    # Convert birth and selected date/time from local time zone to UTC
    birth_date_utc, birth_time_utc = convert_local_to_utc(birth_date, birth_time, latitude, longitude)
    today_date_utc, today_time_utc = convert_local_to_utc(today_ist, fixed_time_ist, custom_latitude, custom_longitude)

    # Get planetary positions
    raw_positions1, formatted_positions1, natal_aspects1 = get_planetary_positions(birth_date_utc, birth_time_utc, latitude, longitude, include_asc_midheaven=True)
    raw_positions2, formatted_positions2, natal_aspects2 = get_planetary_positions(today_date_utc, today_time_utc, custom_latitude, custom_longitude, include_asc_midheaven=True)

    # Generate aspects between birth and comparison planetary positions
    aspects = compare_planetary_positions(raw_positions1, raw_positions2)

    # Generate aspects between birth and comparison planetary positions
    aspects_2 = compare_planetary_positions(raw_positions2, raw_positions1)

    # Generate natal aspects within birth planetary positions
    sorted_natal_aspects = calculate_and_sort_natal_aspects(raw_positions1)

    # Generate natal aspects within comparison planetary positions
    sorted_natal_aspects_2 = calculate_and_sort_natal_aspects(raw_positions2)

    # Calculate house cusps for Person A and Person B
    house_cusps_person_a = calculate_house_cusps(raw_positions1['Ascendant'])
    house_cusps_person_b = calculate_house_cusps(raw_positions2['Ascendant'])

    # Map planets to houses for both persons
    planets_in_houses_a_in_b = map_planets_to_houses(raw_positions1, house_cusps_person_b)
    planets_in_houses_b_in_a = map_planets_to_houses(raw_positions2, house_cusps_person_a)

    # Calculate the Composite Chart
    composite_positions, composite_house_cusps = calculate_composite_chart(raw_positions1, raw_positions2)

    # Calculate aspects within the composite chart
    composite_aspects = calculate_and_sort_natal_aspects(composite_positions)

    # Create an in-memory file
    file_content = io.StringIO()
    file_content.write(f"Planetary positions for Person A ({birth_date} {birth_time} local time):\n")

    for planet, data in formatted_positions1.items():
        house, degrees_into_house = map_planets_to_houses(raw_positions1, house_cusps_person_a).get(planet, (None, None))

        if house is not None:
            file_content.write(f"{planet}: {data} (House {house}, {degrees_into_house:.2f}° into house)\n")
        else:
            file_content.write(f"{planet}: {data}\n")

        # Append aspects for this planet
        if planet in natal_aspects1 and natal_aspects1[planet]:
            file_content.write(f"  Aspects:\n")
            for aspect_info in natal_aspects1[planet]:
                file_content.write(f"    - {aspect_info}\n")


    # Dashas for Person A
    birth_dt_obj = datetime.strptime(birth_date_utc + ' ' + birth_time_utc, '%Y/%m/%d %H:%M:%S')
    dashas_a = get_vimshottari_dashas(raw_positions1['Moon'], birth_dt_obj)

    file_content.write("\nDashas of Person A:\n")
    for maha, start, end, antars in dashas_a:
        file_content.write(f"\n{maha}: {start.strftime('%Y-%m-%d')} to {end.strftime('%Y-%m-%d')}\n")
        for antar, antar_start, antar_end in antars:
            file_content.write(f"  - {antar}: {antar_start.strftime('%Y-%m-%d')} to {antar_end.strftime('%Y-%m-%d')}\n")


    # Prepare Julian date for today and return tracking config
    now_julian = swe.julday(datetime.utcnow().year, datetime.utcnow().month, datetime.utcnow().day)

    planet_table = [
        ("Moon", swe.MOON, 0.05, 13),
        ("Mercury", swe.MERCURY, 0.02, 30),
        ("Venus", swe.VENUS, 0.02, 30),
        ("Sun", swe.SUN, 0.01, 45),
        ("Mars", swe.MARS, 0.02, 20),
        ("Jupiter", swe.JUPITER, 0.02, 15),
        ("Saturn", swe.SATURN, 0.01, 10),
        ("Uranus", swe.URANUS, 0.01, 5),
        ("Neptune", swe.NEPTUNE, 0.01, 5),
        ("Pluto", swe.PLUTO, 0.01, 5),
        ("Rahu", swe.TRUE_NODE, 0.02, 8)
    ]

    # Planetary Return Ages for Person A
    file_content.write("\nPlanetary Return Ages for Person A (as of today):\n")

    year_a, month_a, day_a = map(int, birth_date_utc.split('/'))
    hour_a, minute_a, second_a = map(float, birth_time_utc.split(':'))
    decimal_hour_a = hour_a + minute_a / 60 + second_a / 3600
    birth_julian_a = swe.julday(year_a, month_a, day_a, decimal_hour_a)

    returns_a = {}
    for name, pid, orb, cooldown in planet_table:
        returns_a[name] = calculate_return_age(pid, raw_positions1[name], birth_julian_a, now_julian, orb, cooldown)
        file_content.write(f"{name}: {returns_a[name]:.2f} returns\n")
    file_content.write(f"Ketu: {returns_a['Rahu']:.2f} returns\n")

    file_content.write(f"\nPlanetary positions for Person B ({today_ist} {fixed_time_ist} local time):\n")

    for planet, data in formatted_positions2.items():
        house, degrees_into_house = map_planets_to_houses(raw_positions2, house_cusps_person_b).get(planet, (None, None))

        if house is not None:
            file_content.write(f"{planet}: {data} (House {house}, {degrees_into_house:.2f}° into house)\n")
        else:
            file_content.write(f"{planet}: {data}\n")

        # Append aspects for this planet
        if planet in natal_aspects2 and natal_aspects2[planet]:
            file_content.write(f"  Aspects:\n")
            for aspect_info in natal_aspects2[planet]:
                file_content.write(f"    - {aspect_info}\n")


    # Dashas for Person B
    today_dt_obj = datetime.strptime(today_date_utc + ' ' + today_time_utc, '%Y/%m/%d %H:%M:%S')
    dashas_b = get_vimshottari_dashas(raw_positions2['Moon'], today_dt_obj)

    file_content.write("\nDashas of Person B:\n")
    for maha, start, end, antars in dashas_b:
        file_content.write(f"\n{maha}: {start.strftime('%Y-%m-%d')} to {end.strftime('%Y-%m-%d')}\n")
        for antar, antar_start, antar_end in antars:
            file_content.write(f"  - {antar}: {antar_start.strftime('%Y-%m-%d')} to {antar_end.strftime('%Y-%m-%d')}\n")


    # Planetary Return Ages for Person B
    file_content.write("\nPlanetary Return Ages for Person B (as of today):\n")

    year_b, month_b, day_b = map(int, today_date_utc.split('/'))
    hour_b, minute_b, second_b = map(float, today_time_utc.split(':'))
    decimal_hour_b = hour_b + minute_b / 60 + second_b / 3600
    birth_julian_b = swe.julday(year_b, month_b, day_b, decimal_hour_b)

    returns_b = {}
    for name, pid, orb, cooldown in planet_table:
        returns_b[name] = calculate_return_age(pid, raw_positions2[name], birth_julian_b, now_julian, orb, cooldown)
        file_content.write(f"{name}: {returns_b[name]:.2f} returns\n")
    file_content.write(f"Ketu: {returns_b['Rahu']:.2f} returns\n")

    # Add planet-to-house mappings to the astrology report
    file_content.write("\nPerson A's planets in Person B's houses:\n")

    for planet, (house, degrees_into_house) in planets_in_houses_a_in_b.items():
        file_content.write(f"{planet}: House {house}, {degrees_into_house:.2f}° into house\n")

        # Append synastry aspects involving this planet (Person A to Person B)
        aspect_list = [aspect_info for aspect_info in aspects if aspect_info[0] == planet]

        if aspect_list:
            file_content.write(f"  Synastry Aspects:\n")
            for aspect_info in aspect_list:
                file_content.write(f"    - {aspect_info[2]} {aspect_info[1]} (Orb: {aspect_info[3]:.2f}°)\n")

    file_content.write("\nPerson B's planets in Person A's houses:\n")

    for planet, (house, degrees_into_house) in planets_in_houses_b_in_a.items():
        file_content.write(f"{planet}: House {house}, {degrees_into_house:.2f}° into house\n")

        # Append synastry aspects involving this planet (Person B to Person A)
        aspect_list = [aspect_info for aspect_info in aspects_2 if aspect_info[0] == planet]

        if aspect_list:
            file_content.write(f"  Synastry Aspects:\n")
            for aspect_info in aspect_list:
                file_content.write(f"    - {aspect_info[2]} {aspect_info[1]} (Orb: {aspect_info[3]:.2f}°)\n")

    file_content.write("\nComposite Chart of Person A and Person B:\n")

    for planet, degree in composite_positions.items():
        house, degrees_into_house = map_planets_to_houses(composite_positions, composite_house_cusps).get(planet, (None, None))
        zodiac_info = f"{get_zodiac(degree)}, {get_nakshatra(degree)}"

        # Write planet's position with house placement
        if house is not None:
            file_content.write(f"{planet}: {degree:.2f}° ({zodiac_info}) (House {house}, {degrees_into_house:.2f}° into house)\n")
        else:
            file_content.write(f"{planet}: {degree:.2f}° ({zodiac_info})\n")

        # Append aspects for this planet
        if planet in composite_aspects and composite_aspects[planet]:
            file_content.write(f"  Aspects:\n")
            for aspect_info in composite_aspects[planet]:
                file_content.write(f"    - {aspect_info}\n")

    file_content.write("\nHow to Use Planetary Keywords to Interpret a Person\n")

    keywords = {
        "Ascendant": ("How the person appears and operates in the world",
                    "Represents the person’s outward behavior, physical appearance, and default life approach."),
        "Ketu": ("Where the person is instinctively wise but detached",
                "Indicates areas of unconscious intelligence, past-life mastery, and non-attachment."),
        "Venus": ("How the person loves, attracts, and finds beauty",
                "Reveals how they engage in relationships, what they value in love and aesthetics, and how they experience pleasure."),
        "Sun": ("The person’s core identity and life purpose",
                "Represents their ego, confidence, and conscious self."),
        "Moon": ("The person’s emotional nature and inner life",
                "Shows their emotional needs, memory, and subconscious responses."),
        "Mars": ("How the person acts, desires, and asserts themselves",
                "Governs energy, drive, sexuality, and confidence."),
        "Rahu": ("The person’s hunger, obsession, and ambition",
                "Points to desires in this lifetime — where they’re obsessed, restless, or seeking material validation."),
        "Jupiter": ("The person’s faith, expansion, and wisdom",
                    "Reveals areas of natural growth, optimism, and teaching."),
        "Saturn": ("Where the person must grow through effort and structure",
                "Describes their relationship to discipline, fear, and responsibility."),
        "Mercury": ("How the person thinks, learns, and communicates",
                    "Covers mental habits, intellect, and expression.")
    }

    for planet, (theme, description) in keywords.items():
        file_content.write(f"\n{planet} – \"{theme}\"\n{description}\n")

    file_content.write("\nHow to Interpret Aspects – Quick Guide\n")

    aspects_guide = {
        "Conjunction (0°)": "Merge the energies - The planets fuse into one strong expression. Look for intensity, unity, or potential overload. Use them as a single force.",
        "Opposition (180°)": "Balance the energies - Polarity needs integration. Don’t choose sides — work toward awareness, perspective, and cooperation.",
        "Trine (120°)": "Trust and apply the ease - Natural talent or harmony. Don’t coast — respect and use the flow for growth or mastery.",
        "Square (90°)": "Adjust through action - Tension needs movement. This is where work is required. Change, break patterns, build resilience.",
        "Sextile (60°)": "Activate the opportunity - Potential is there, but it needs a choice or effort to become real. Take initiative.",
        "Quincunx (150°)": "Realign what's out of sync - Mismatch creates discomfort. Something must shift — change habits, health, or perspective.",
        "Semi-Square (45°)": "Notice subtle friction - Irritation that builds. Pay attention to low-level stress and act before it escalates.",
        "Sesquiquadrate (135°)": "Work through background pressure - Long-term tension beneath the surface. Be mindful of patterns that quietly sabotage.",
        "Semi-Sextile (30°)": "Gently adjust and adapt - Not obvious, but valuable. Shows awkward alliances or emerging insight — lean into the nuance."
    }

    for aspect, meaning in aspects_guide.items():
        file_content.write(f"\n{aspect}\n{meaning}\n")


    file_content.write("\nHow to Use House Keywords to Interpret a Person\n")

    houses_guide = {
        "1st House": ("How the person shows up in life", "Personality, approach, physical body, overall life direction."),
        "2nd House": ("What the person values and builds", "Finances, possessions, talents, self-worth."),
        "3rd House": ("How the person thinks, speaks, and connects locally", "Mindset, communication, siblings, short travel."),
        "4th House": ("Where the person comes from emotionally and ancestrally", "Family, home, roots, emotional security."),
        "5th House": ("How the person creates, plays, and seeks joy", "Romance, creativity, children, hobbies."),
        "6th House": ("How the person manages routine and responsibility", "Work, health, service, daily habits."),
        "7th House": ("How the person partners and reflects through others", "Marriage, contracts, close relationships, projection."),
        "8th House": ("Where the person transforms, merges, and releases control", "Death, sex, shared resources, deep psychology."),
        "9th House": ("How the person expands their world and beliefs", "Travel, higher learning, spirituality, philosophy."),
        "10th House": ("What the person aims to achieve and be known for", "Career, public life, reputation, authority."),
        "11th House": ("How the person aligns with groups and visions", "Friendships, community, future goals, ideals."),
        "12th House": ("Where the person withdraws, hides, or transcends", "Isolation, dreams, karma, spiritual retreat, the unconscious."),
    }

    for house, (theme, description) in houses_guide.items():
        file_content.write(f"\n{house} – \"{theme}\"\n{description}\n")

    file_content.write("\nHow to Interpret Planetary Lordships\n")
    file_content.write("Rule:\nPlanet rules one house, acts in another.\n")
    file_content.write("The ruler’s house shows where that house’s story unfolds.\n")
    file_content.write("The sign it’s in shows how it behaves.\n")
    file_content.write("\nInterpret like this:\n")
    file_content.write("• \"Lord of the Xth in the Yth\" → X house issues play out in Y house matters.\n")
    file_content.write("\nExamples:\n")
    file_content.write("• Lord of the 1st in the 4th → Self expresses through home and emotional roots.\n")
    file_content.write("• Lord of the 10th in the 6th → Career plays out through work ethic, service, or health field.\n")
    file_content.write("• Lord of the 7th in the 2nd → Relationships affect income and values.\n")
    file_content.write("\nKeep in mind:\n")
    file_content.write("• House being ruled = theme\n")
    file_content.write("• House occupied = stage\n")
    file_content.write("• Sign = style or tone\n")

    print("Astrology report created in memory.")


if __name__ == "__main__":
    main()


def run_astrology_report(birth_date, birth_time, latitude, longitude,
                         compare_date, compare_time,
                         compare_latitude, compare_longitude):
    import io
    from datetime import datetime
    import swisseph as swe

    # Convert to UTC
    birth_date_utc, birth_time_utc = convert_local_to_utc(birth_date, birth_time, latitude, longitude)
    compare_date_utc, compare_time_utc = convert_local_to_utc(compare_date, compare_time, compare_latitude, compare_longitude)

    # Get planetary data
    raw_positions1, formatted_positions1, natal_aspects1 = get_planetary_positions(birth_date_utc, birth_time_utc, latitude, longitude, include_asc_midheaven=True)
    raw_positions2, formatted_positions2, natal_aspects2 = get_planetary_positions(compare_date_utc, compare_time_utc, compare_latitude, compare_longitude, include_asc_midheaven=True)
    aspects = compare_planetary_positions(raw_positions1, raw_positions2)
    aspects_2 = compare_planetary_positions(raw_positions2, raw_positions1)

    house_cusps_a = calculate_house_cusps(raw_positions1['Ascendant'])
    house_cusps_b = calculate_house_cusps(raw_positions2['Ascendant'])
    in_houses_a_in_b = map_planets_to_houses(raw_positions1, house_cusps_b)
    in_houses_b_in_a = map_planets_to_houses(raw_positions2, house_cusps_a)

    composite_positions, composite_house_cusps = calculate_composite_chart(raw_positions1, raw_positions2)
    composite_aspects = calculate_and_sort_natal_aspects(composite_positions)

    # === Section 1: Birth Chart A ===
    chart_a = io.StringIO()
    chart_a.write(f"Planetary positions for Person A ({birth_date} {birth_time}):\n")
    for planet, data in formatted_positions1.items():
        house, deg = map_planets_to_houses(raw_positions1, house_cusps_a).get(planet, (None, None))
        if house:
            chart_a.write(f"{planet}: {data} (House {house}, {deg:.2f}° into house)\n")
        else:
            chart_a.write(f"{planet}: {data}\n")
        for aspect in natal_aspects1.get(planet, []):
            chart_a.write(f"  Aspect: {aspect}\n")

    # Dashas A
    birth_dt_obj = datetime.strptime(f"{birth_date_utc} {birth_time_utc}", "%Y/%m/%d %H:%M:%S")
    dashas_a = get_vimshottari_dashas(raw_positions1['Moon'], birth_dt_obj)

    # Return Ages A
    now_julian = swe.julday(datetime.utcnow().year, datetime.utcnow().month, datetime.utcnow().day)
    y, m, d = map(int, birth_date_utc.split('/'))
    h, mi, s = map(float, birth_time_utc.split(':'))
    decimal_hour = h + mi / 60 + s / 3600
    birth_julian = swe.julday(y, m, d, decimal_hour)

    returns_a = io.StringIO()
    returns_a.write("Dashas for Person A:\n")
    for maha, start, end, antars in dashas_a:
        returns_a.write(f"\n{maha}: {start.date()} to {end.date()}\n")
        for antar, s, e in antars:
            returns_a.write(f"  - {antar}: {s.date()} to {e.date()}\n")

    returns_a.write("\nPlanetary Return Ages:\n")
    planet_table = [("Moon", swe.MOON, 0.05, 13), ("Mercury", swe.MERCURY, 0.02, 30),
                    ("Venus", swe.VENUS, 0.02, 30), ("Sun", swe.SUN, 0.01, 45),
                    ("Mars", swe.MARS, 0.02, 20), ("Jupiter", swe.JUPITER, 0.02, 15),
                    ("Saturn", swe.SATURN, 0.01, 10), ("Uranus", swe.URANUS, 0.01, 5),
                    ("Neptune", swe.NEPTUNE, 0.01, 5), ("Pluto", swe.PLUTO, 0.01, 5),
                    ("Rahu", swe.TRUE_NODE, 0.02, 8)]
    for name, pid, orb, cooldown in planet_table:
        ret = calculate_return_age(pid, raw_positions1[name], birth_julian, now_julian, orb, cooldown)
        returns_a.write(f"{name}: {ret:.2f} returns\n")
    returns_a.write(f"Ketu: {ret:.2f} returns\n")

    # === Section 2: Birth Chart B ===
    chart_b = io.StringIO()
    chart_b.write(f"Planetary positions for Person B ({compare_date} {compare_time}):\n")
    for planet, data in formatted_positions2.items():
        house, deg = map_planets_to_houses(raw_positions2, house_cusps_b).get(planet, (None, None))
        if house:
            chart_b.write(f"{planet}: {data} (House {house}, {deg:.2f}° into house)\n")
        else:
            chart_b.write(f"{planet}: {data}\n")
        for aspect in natal_aspects2.get(planet, []):
            chart_b.write(f"  Aspect: {aspect}\n")

    # Dashas and Returns B
    today_dt_obj = datetime.strptime(f"{compare_date_utc} {compare_time_utc}", "%Y/%m/%d %H:%M:%S")
    dashas_b = get_vimshottari_dashas(raw_positions2['Moon'], today_dt_obj)
    yb, mb, db = map(int, compare_date_utc.split('/'))
    hb, mib, sb = map(float, compare_time_utc.split(':'))
    birth_julian_b = swe.julday(yb, mb, db, hb + mib / 60 + sb / 3600)

    returns_b = io.StringIO()
    returns_b.write("Dashas for Person B:\n")
    for maha, start, end, antars in dashas_b:
        returns_b.write(f"\n{maha}: {start.date()} to {end.date()}\n")
        for antar, s, e in antars:
            returns_b.write(f"  - {antar}: {s.date()} to {e.date()}\n")

    returns_b.write("\nPlanetary Return Ages:\n")
    for name, pid, orb, cooldown in planet_table:
        ret = calculate_return_age(pid, raw_positions2[name], birth_julian_b, now_julian, orb, cooldown)
        returns_b.write(f"{name}: {ret:.2f} returns\n")
    returns_b.write(f"Ketu: {ret:.2f} returns\n")

    # === Section 3: Composite Chart ===
    composite = io.StringIO()
    composite.write("Composite Chart:\n")
    for planet, degree in composite_positions.items():
        house, into = map_planets_to_houses(composite_positions, composite_house_cusps).get(planet, (None, None))
        zodiac_info = f"{get_zodiac(degree)}, {get_nakshatra(degree)}"
        if house:
            composite.write(f"{planet}: {degree:.2f}° ({zodiac_info}) (House {house}, {into:.2f}° into house)\n")
        else:
            composite.write(f"{planet}: {degree:.2f}° ({zodiac_info})\n")
        for asp in composite_aspects.get(planet, []):
            composite.write(f"  Aspect: {asp}\n")

    # === Section 4: Synastry ===
    synastry = io.StringIO()
    synastry.write("Person A's planets in Person B's houses:\n")
    for planet, (house, deg) in in_houses_a_in_b.items():
        synastry.write(f"{planet}: House {house}, {deg:.2f}° into house\n")
        for asp in [a for a in aspects if a[0] == planet]:
            synastry.write(f"  Synastry Aspect: {asp[2]} {asp[1]} (Orb: {asp[3]:.2f}°)\n")

    synastry.write("\nPerson B's planets in Person A's houses:\n")
    for planet, (house, deg) in in_houses_b_in_a.items():
        synastry.write(f"{planet}: House {house}, {deg:.2f}° into house\n")
        for asp in [a for a in aspects_2 if a[0] == planet]:
            synastry.write(f"  Synastry Aspect: {asp[2]} {asp[1]} (Orb: {asp[3]:.2f}°)\n")

    # === Return all sections as structured output ===
    return {
        "birth_chart_a": chart_a.getvalue(),
        "birth_chart_b": chart_b.getvalue(),
        "composite_summary": composite.getvalue(),
        "dashas": returns_a.getvalue() + "\n\n" + returns_b.getvalue(),
        "compatibility": synastry.getvalue()
    }