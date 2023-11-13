import pandas as pd
import matplotlib.pyplot as plt
import os
import plotly.graph_objects as go
import json
import plotly

FOLDER_PATH = "static/plots_and_csvs"


def indicate_by_threshold(x, threshold):
    return 1 if x > threshold else 0


def plot_custom_plot(df, current_time, threshold, words, index=0):
    """Get the output dataframe from the model and manipulate in to a plot"""
    try:
        df.columns = ["protein", "index", "amino acid", "score", "default indicator", "threshold"]
    except:
        df.columns = ["protein", "index", "amino acid", "score", "default indicator", "threshold", "custom threshold", "custom indicator"]

    df = df[df["protein"].astype(str) == str(words[index])]

    if not threshold:
        threshold = df.threshold.iloc[0]
    else:
        threshold = float(threshold)

    x_axis = [i for i in range(df.shape[0])]
    y_axis = list(df.score.astype(float))
    custom_labels = [df["amino acid"].iloc[i] for i in x_axis]
    custom_labels = [custom_labels[i].upper() if y_axis[i] > threshold else custom_labels[i].lower() for i in
                     range(len(custom_labels))]

    x_bad = [i for i in range(len(y_axis)) if y_axis[i] > threshold]
    y_bad = [y for y in y_axis if y > threshold]
    custom_bad = [custom_labels[i] for i in x_bad]

    trace1 = go.Scatter(x=x_axis, y=y_axis, mode='lines',
                        name='scores', hovertext=custom_labels, line=dict(color="turquoise"), showlegend=True)
    trace2 = go.Scatter(x=x_axis, y=[threshold for i in range(df.shape[0])],
                        mode='lines', name='threshold', line=dict(color="blue", dash="dot"), showlegend=True)
    trace3 = go.Scatter(x=x_bad, y=y_bad, mode='markers', name="epitope", hovertext=custom_bad,
                        line=dict(color="salmon"), showlegend=True)

    fig = go.Figure(data=[trace2, trace1, trace3])

    fig.update_layout(
        title=f"{words[index]}",
        title_x=0.5,
        xaxis_title="Residue Position",
        yaxis_title="Score",
        legend_title="Legend",
        legend_x=1.0,  # Adjust the x position
        legend_y=1.0,  # Adjust the y position
        # legend_xanchor="center",  # Set the x anchor point
        legend_bgcolor="lightgray",  # Set background color
        legend_bordercolor="gray",  # Set border color
    )

    save_path = os.path.join(FOLDER_PATH, "scores_plot_" + str(current_time) + ".json")
    plotly_data = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    # with open(save_path, "w") as f:
    #     json.dump(plotly_data, f)
    # fig.show()

    return plotly_data, save_path, threshold


def create_new_csv(df, current_time, threshold):
    """Create a new custom csv file"""
    # original_df_path = os.path.join("plots_and_csvs", "model_output"+str(current_time)+".csv")
    original_df_path = os.path.join(FOLDER_PATH, "model_output" + str(current_time) + ".csv")
    try:
        df.columns = ["protein", "index", "amino acid", "score", "default indicator", "default threshold"]
    except:
        df.columns = ["protein", "index", "amino acid", "score", "default indicator", "default threshold", "custom threshold", "custom indicator"]
    if threshold:
        threshold = float(threshold)
        df["custom threshold"] = threshold
        df["custom indicator"] = df["score"].apply(indicate_by_threshold, threshold=threshold)
    df.to_csv(original_df_path, index=False)

    proteins = list(set(list(df["protein"])))
    words = []

    if threshold:
        for protein in proteins:
            new_word = ""
            restricted_df = df[df["protein"] == protein]
            for i in range(restricted_df.shape[0]):
                acid = restricted_df["amino acid"].iloc[i]
                new_word += acid.upper() if restricted_df["custom indicator"].iloc[i] == 1 else acid.lower()
            words.append(new_word)
    else:
        for protein in proteins:
            new_word = ""
            restricted_df = df[df["protein"] == protein]
            for i in range(restricted_df.shape[0]):
                acid = restricted_df["amino acid"].iloc[i]
                new_word += acid.upper() if restricted_df["default indicator"].iloc[i] == 1 else acid.lower()
            words.append(new_word)

    # new_df = pd.DataFrame({"protein": proteins, "sequences": words})
    # # save_path = os.path.join("plots_and_csvs", "protein_epitope_upper_case"+str(current_time)+".csv")
    save_path = os.path.join(FOLDER_PATH, "protein_epitope_upper_case" + str(current_time) + ".fasta")
    # new_df.to_csv(save_path)

    with open(save_path, 'w') as f:
        for i in range(len(proteins)):
            f.write(">"+proteins[i])
            f.write("\n")
            f.write(words[i])
            f.write("\n")

    return save_path, original_df_path, proteins
