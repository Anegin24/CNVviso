from flask import Flask, render_template, request
from tkinter import filedialog
import subprocess
import os

app = Flask(__name__)

# Function to execute the Bash script with the selected directory
def execute_script(input_directory):
    if not os.path.isdir(input_directory):
        return "Invalid directory path"

    # Get the directory of the Python script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    
    # Construct the path to the Bash script relative to the Python script
    bash_script_path = os.path.join(script_dir, "test.sh")

    # Capture the output of the Bash script
    process = subprocess.Popen(['bash', bash_script_path, input_directory], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    
    # Combine stdout and stderr
    output = stdout.decode('utf-8') + '\n' + stderr.decode('utf-8')

    # Return the output
    return output

@app.route('/', methods=['GET', 'POST'])
def index():
    output = None
    if request.method == 'POST':
        directory_path = filedialog.askdirectory()
        output = execute_script(directory_path)
    return render_template('index.html', output=output)

if __name__ == '__main__':
    app.run(debug=True)
