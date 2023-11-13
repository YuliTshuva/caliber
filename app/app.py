import datetime
import random
from zipfile import ZipFile
import os
from flask import Flask, request, render_template, flash
import traceback
from run_model import run_model_1, run_model_2, run_model_3
import pandas as pd
import time
from imputations import plot_custom_plot, create_new_csv
import uuid

FOLDER_NAME = "epitope_b_cells_predictor"

random_hex = uuid.uuid4().hex


def get_time():
    """return the current time in seconds"""
    return time.time()


app = Flask(__name__, static_folder='static', template_folder="templates")

app.config['UPLOAD_FOLDER'] = os.path.abspath("upload_folder")

app.config["PERMANENT_SESSION_LIFETIME"] = datetime.timedelta(hours=1)

app.config["SECRET_KEY"] = random_hex


@app.route('/impute-form', methods=['POST'])
def impute_form():
    try:
        # get the input data from the form
        option = request.form.get("option")

        # get the input data from the form
        model = request.form.get("model")

        # get the input data from the form
        encoding = request.form.get("encoding")

        # get the input data from the form
        epitope = request.form.get("epitope")

        # get the input data from the form
        threshold = request.form.get("threshold", "")

        # get the input data from the form
        input_string = request.form.get('seq', "")

        # get the uploaded file from the request object
        file = request.files['attached-file']

        try:
            float(threshold)
        except:
            if threshold != "":
                traceback.print_exc()
                return render_template("error.html", active="", error="The threshold is not numeric.")

        if option == "Protein Sequences" and not input_string and not file:
            return render_template("error.html", active="", error="No protein sequences / PDB ID was entered.")

        # Setting a path to save file
        current_time = get_time()
        save_path = os.path.join(FOLDER_NAME, "input", str(current_time)+".txt")
        extract_path = os.path.join(FOLDER_NAME, "input", str(current_time))
        save_zip_path = os.path.join(FOLDER_NAME, "input", str(current_time)+".zip")

        if option == "Protein Sequences":
            if file:
                file.save(save_path)
            else:
                input_string = input_string.replace(" ", "\n")
                with open(save_path, "w") as f:
                    f.write(input_string)
            run_model_1(model, encoding, epitope, save_path)

        if (option == "PDB IDs" or option == "PDB file") and not file:
            return render_template("error.html", active="", error="No file was uploaded.")

        if option == "PDB IDs":
            # save the file
            file.save(save_path)
            run_model_2(model, encoding, epitope, save_path)

        if option == "PDB file":
            # Save the zip file
            file.save(save_zip_path)
            # Extract the file content
            with ZipFile(save_zip_path, 'r') as zip_file:
                zip_file.extractall(extract_path)
            run_model_3(model, encoding, epitope, extract_path)

        results_path = f"{FOLDER_NAME}/output/predict_predict_{model}_{encoding}_None_{epitope}.csv"
        results_df = pd.read_csv(results_path)

        new_df_path, original_df_path, proteins = create_new_csv(results_df, current_time, threshold)

        plotly_data, plot_path, threshold = plot_custom_plot(results_df, current_time, threshold, proteins)

        # render a template with the imputation results
        return render_template("results_second_edition.html", active="results", threshold=threshold,
                               current_time=current_time, words=proteins, words_send=" ".join(proteins), index=0,
                               plotly_data=plotly_data, custom_df_path=new_df_path, original_df_path=original_df_path,
                               certain_word=proteins[0])

    #render an error template if an exception occurs
    except Exception as e:
        traceback.print_exc()
        return render_template("error.html", active="", error=str(e))


@app.route("/update-threshold", methods=["POST"])
def update_threshold():
    # get the input data from the form
    threshold = request.form.get("threshold", "")

    try:
        float(threshold)
    except:
        traceback.print_exc()
        return render_template("error.html", active="", error="The threshold is not numeric.")

    threshold = float(threshold)
    if threshold > 0.5:
        threshold = 0.5
    if threshold < 0:
        threshold = 0

    original_df_path = request.form.get("original_df_path")
    current_time = request.form.get("current_time")
    words_send = request.form.get("words_send2")
    words = words_send.split()
    index = int(request.form.get("index"))

    results_df = pd.read_csv(original_df_path)

    plotly_data, _, _ = plot_custom_plot(results_df, current_time, threshold, words)

    new_df_path, original_df_path, _ = create_new_csv(results_df, current_time, threshold)

    return render_template("results_second_edition.html", active="results", threshold=threshold, words=words,
                           words_send=words_send, plotly_data=plotly_data, custom_df_path=new_df_path, index=index,
                           original_df_path=original_df_path, current_time=current_time, certain_word=words[index])


@app.route('/update-plot', methods=['POST'])
def update_plot():
    protein = request.form.get("protein-val")

    words_send = request.form.get("words_send")
    words = words_send.split()

    original_df_path = request.form.get("df_path")
    current_time = request.form.get("current_time_2")
    threshold = request.form.get("threshold-val")
    new_df_path = request.form.get("new_df_path")
    index = words.index(protein)

    plotly_data, _, _ = plot_custom_plot(pd.read_csv(original_df_path),
                                         current_time,
                                         threshold,
                                         words,
                                         index)

    return render_template("results_second_edition.html", active="results", threshold=threshold, words=words,
                           words_send=words_send, plotly_data=plotly_data, custom_df_path=new_df_path, index=index,
                           original_df_path=original_df_path, current_time=current_time, certain_word=words[index])


@app.route('/', methods=['GET'])
@app.route('/Home', methods=['GET'])
def home():
    return render_template("index.html", active="Home")


@app.route('/Example', methods=['GET'])
def example():
    return render_template("example.html", active="Example")


@app.route('/Help', methods=['GET'])
def help():
    return render_template("help.html", active="Help")


@app.route('/About', methods=['GET'])
def about():
    return render_template("about.html", active="About")
#
#
# @app.route('/experiments', methods=['GET', 'POST'])
# def experiments():
#     with open("static\plots_and_csvs\scores_plot_1698302896.858246.json", "r") as f:
#         data = json.load(f)
#     return render_template("results_second_edition.html", active="results",
#                            plotly_data=data,
#                            custom_df_path="static\plots_and_csvs\protein_epitope_upper_case1698302896.858246.csv",
#                            original_df_path="static\plots_and_csvs\model_output1698302896.858246.csv")


if __name__ == "__main__":
    app.run(debug=True, host='0.0.0.0', port=5000, use_reloader=True)
