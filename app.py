from flask import Flask, render_template, request
# Placeholder for your astrology logic (you'll wire this in later)
import main

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def home():
    result = None
    if request.method == 'POST':
        # Get form input values
        birth_date = request.form['birth_date']
        birth_time = request.form['birth_time']
        latitude = request.form['latitude']
        longitude = request.form['longitude']
        compare_date = request.form['compare_date']
        compare_time = request.form['compare_time']
        compare_latitude = request.form['compare_latitude']
        compare_longitude = request.form['compare_longitude']

        # For now, just combine the input to test if it's working
        result = main.run_astrology_report(
            birth_date, birth_time,
            latitude, longitude,
            compare_date, compare_time,
            compare_latitude, compare_longitude
        )

        # Later we'll call main.run_astrology_report(...) instead of just combining text

    return render_template('form.html', result_sections=result)

if __name__ == '__main__':
    app.run(debug=True)