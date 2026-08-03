import marimo

__generated_with = "0.23.9"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _(mo):
    mo.md(
        "# Transfection CSV plots\n\n"
        "Plain pandas + matplotlib. Edit the CSV paths below if needed."
    )
    return


@app.cell
def _(mo):
    import pandas as pd

    sc0 = pd.read_csv(r"C:\Users\ctyja\data\20260603\timeseries\sc0_ch1.csv")
    sc1 = pd.read_csv(r"C:\Users\ctyja\data\20260603\timeseries\sc1_ch1.csv")
    sc2 = pd.read_csv(r"C:\Users\ctyja\data\20260603\timeseries\sc2_ch1.csv")
    sc3 = pd.read_csv(r"C:\Users\ctyja\data\20260603\timeseries\sc3_ch1.csv")
    sc4 = pd.read_csv(r"C:\Users\ctyja\data\20260603\timeseries\sc4_ch1.csv")
    sc5 = pd.read_csv(r"C:\Users\ctyja\data\20260603\timeseries\sc5_ch1.csv")
    auc_df = pd.read_csv(r"C:\Users\ctyja\data\20260603\results\auc.csv")
    fit_df = pd.read_csv(r"C:\Users\ctyja\data\20260603\results\fit.csv")

    mo.vstack(
        [
            mo.md("### sc0_ch1.csv"),
            mo.ui.table(sc0.head(10)),
            mo.md("### auc.csv"),
            mo.ui.table(auc_df),
            mo.md("### fit.csv"),
            mo.ui.table(fit_df),
        ]
    )
    return auc_df, fit_df, sc0, sc1, sc2, sc3, sc4, sc5


@app.cell
def _(auc_df, fit_df, mo, sc0, sc1, sc2, sc3, sc4, sc5):
    import matplotlib.pyplot as plt
    import numpy as np

    fig, axes = plt.subplots(2, 3, figsize=(12.0, 8.0), squeeze=False)

    ax = axes[0, 0]
    for _, trace in sc0.groupby(["pos", "roi"], sort=True):
        ax.plot(
            trace["t"].to_numpy(dtype=float) * 10.0,
            trace["corrected"],
            color="green",
            alpha=0.1,
        )
    ax.set_ylim(
        np.percentile(sc0["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc0["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc0_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[0, 1]
    for _, trace in sc1.groupby(["pos", "roi"], sort=True):
        ax.plot(
            trace["t"].to_numpy(dtype=float) * 10.0,
            trace["corrected"],
            color="green",
            alpha=0.1,
        )
    ax.set_ylim(
        np.percentile(sc1["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc1["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc1_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[0, 2]
    for _, trace in sc2.groupby(["pos", "roi"], sort=True):
        ax.plot(
            trace["t"].to_numpy(dtype=float) * 10.0,
            trace["corrected"],
            color="green",
            alpha=0.1,
        )
    ax.set_ylim(
        np.percentile(sc2["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc2["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc2_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[1, 0]
    for _, trace in sc3.groupby(["pos", "roi"], sort=True):
        ax.plot(
            trace["t"].to_numpy(dtype=float) * 10.0,
            trace["corrected"],
            color="green",
            alpha=0.1,
        )
    ax.set_ylim(
        np.percentile(sc3["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc3["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc3_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[1, 1]
    for _, trace in sc4.groupby(["pos", "roi"], sort=True):
        ax.plot(
            trace["t"].to_numpy(dtype=float) * 10.0,
            trace["corrected"],
            color="green",
            alpha=0.1,
        )
    ax.set_ylim(
        np.percentile(sc4["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc4["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc4_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[1, 2]
    for _, trace in sc5.groupby(["pos", "roi"], sort=True):
        ax.plot(
            trace["t"].to_numpy(dtype=float) * 10.0,
            trace["corrected"],
            color="green",
            alpha=0.1,
        )
    ax.set_ylim(
        np.percentile(sc5["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc5["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc5_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    fig.tight_layout()
    traces_fig = fig

    fig, ax = plt.subplots(figsize=(12.0, 8.0))
    ax.boxplot(
        [
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 0), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 1), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 2), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 3), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 4), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 5), "auc"
            ].to_numpy(dtype=float),
        ],
        tick_labels=["0", "1", "2", "3", "4", "5"],
    )
    ax.set_xlabel("slide channel")
    ax.set_ylabel("AUC")
    ax.set_ylim(
        np.percentile(
            auc_df.loc[auc_df["auc"] > 0, "auc"].to_numpy(dtype=float), 5
        ),
        np.percentile(
            auc_df.loc[auc_df["auc"] > 0, "auc"].to_numpy(dtype=float), 95
        ),
    )
    fig.tight_layout()
    auc_linear_fig = fig

    fig, ax = plt.subplots(figsize=(12.0, 8.0))
    ax.boxplot(
        [
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 0), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 1), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 2), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 3), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 4), "auc"
            ].to_numpy(dtype=float),
            auc_df.loc[
                (auc_df["auc"] > 0) & (auc_df["slide_channel"] == 5), "auc"
            ].to_numpy(dtype=float),
        ],
        tick_labels=["0", "1", "2", "3", "4", "5"],
    )
    ax.set_xlabel("slide channel")
    ax.set_ylabel("AUC")
    ax.set_yscale("log")
    fig.tight_layout()
    auc_log_fig = fig

    fig, ax = plt.subplots(figsize=(12.0, 8.0))
    ax.boxplot(
        [
            (1.0 / fit_df.loc[
                fit_df["success"].astype(str).str.lower().eq("true")
                & (fit_df["slide_channel"] == 0),
                "protein_decay_rate",
            ]).to_numpy(dtype=float),
            (1.0 / fit_df.loc[
                fit_df["success"].astype(str).str.lower().eq("true")
                & (fit_df["slide_channel"] == 1),
                "protein_decay_rate",
            ]).to_numpy(dtype=float),
            (1.0 / fit_df.loc[
                fit_df["success"].astype(str).str.lower().eq("true")
                & (fit_df["slide_channel"] == 2),
                "protein_decay_rate",
            ]).to_numpy(dtype=float),
            (1.0 / fit_df.loc[
                fit_df["success"].astype(str).str.lower().eq("true")
                & (fit_df["slide_channel"] == 3),
                "protein_decay_rate",
            ]).to_numpy(dtype=float),
            (1.0 / fit_df.loc[
                fit_df["success"].astype(str).str.lower().eq("true")
                & (fit_df["slide_channel"] == 4),
                "protein_decay_rate",
            ]).to_numpy(dtype=float),
            (1.0 / fit_df.loc[
                fit_df["success"].astype(str).str.lower().eq("true")
                & (fit_df["slide_channel"] == 5),
                "protein_decay_rate",
            ]).to_numpy(dtype=float),
        ],
        tick_labels=["0", "1", "2", "3", "4", "5"],
    )
    ax.set_xlabel("slide channel")
    ax.set_ylabel("protein lifetime")
    fig.tight_layout()
    protein_lifetime_fig = fig

    fig, ax = plt.subplots(figsize=(12.0, 8.0))
    ax.boxplot(
        [
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 0),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 0),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 0),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 1),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 1),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 1),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 2),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 2),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 2),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 3),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 3),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 3),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 4),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 4),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 4),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 5),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 5),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 5),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
        ],
        tick_labels=["0", "1", "2", "3", "4", "5"],
    )
    ax.set_xlabel("slide channel")
    ax.set_ylabel("expression rate")
    fig.tight_layout()
    expression_rate_linear_fig = fig

    fig, ax = plt.subplots(figsize=(12.0, 8.0))
    ax.boxplot(
        [
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 0),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 0),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 0),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 1),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 1),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 1),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 2),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 2),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 2),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 3),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 3),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 3),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 4),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 4),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 4),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
            (
                fit_df.loc[
                    fit_df["success"].astype(str).str.lower().eq("true")
                    & (fit_df["slide_channel"] == 5),
                    "expression_amplitude",
                ].to_numpy(dtype=float)
                * (
                    fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 5),
                        "mrna_decay_rate",
                    ].to_numpy(dtype=float)
                    - fit_df.loc[
                        fit_df["success"].astype(str).str.lower().eq("true")
                        & (fit_df["slide_channel"] == 5),
                        "protein_decay_rate",
                    ].to_numpy(dtype=float)
                )
            ),
        ],
        tick_labels=["0", "1", "2", "3", "4", "5"],
    )
    ax.set_xlabel("slide channel")
    ax.set_ylabel("expression rate")
    ax.set_yscale("log")
    fig.tight_layout()
    expression_rate_log_fig = fig

    fig, axes = plt.subplots(2, 3, figsize=(12.0, 8.0), squeeze=False)

    ax = axes[0, 0]
    for (pos, roi), trace in sc0.groupby(["pos", "roi"], sort=True):
        row = fit_df.loc[
            (fit_df["success"].astype(str).str.lower().eq("true"))
            & (fit_df["slide_channel"] == 0)
            & (fit_df["pos"] == pos)
            & (fit_df["roi"] == roi)
        ]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        minutes = trace["t"].to_numpy(dtype=float) * 10.0
        dt = np.maximum(minutes - row["onset_time"], 0.0)
        predicted = row["baseline_intensity"] + row["expression_amplitude"] * (
            np.exp(-row["protein_decay_rate"] * dt)
            - np.exp(-row["mrna_decay_rate"] * dt)
        )
        predicted = np.where(
            minutes < row["onset_time"], row["baseline_intensity"], predicted
        )
        ax.plot(minutes, predicted, color="green", alpha=0.1)
    ax.set_ylim(
        np.percentile(sc0["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc0["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc0_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[0, 1]
    for (pos, roi), trace in sc1.groupby(["pos", "roi"], sort=True):
        row = fit_df.loc[
            (fit_df["success"].astype(str).str.lower().eq("true"))
            & (fit_df["slide_channel"] == 1)
            & (fit_df["pos"] == pos)
            & (fit_df["roi"] == roi)
        ]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        minutes = trace["t"].to_numpy(dtype=float) * 10.0
        dt = np.maximum(minutes - row["onset_time"], 0.0)
        predicted = row["baseline_intensity"] + row["expression_amplitude"] * (
            np.exp(-row["protein_decay_rate"] * dt)
            - np.exp(-row["mrna_decay_rate"] * dt)
        )
        predicted = np.where(
            minutes < row["onset_time"], row["baseline_intensity"], predicted
        )
        ax.plot(minutes, predicted, color="green", alpha=0.1)
    ax.set_ylim(
        np.percentile(sc1["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc1["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc1_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[0, 2]
    for (pos, roi), trace in sc2.groupby(["pos", "roi"], sort=True):
        row = fit_df.loc[
            (fit_df["success"].astype(str).str.lower().eq("true"))
            & (fit_df["slide_channel"] == 2)
            & (fit_df["pos"] == pos)
            & (fit_df["roi"] == roi)
        ]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        minutes = trace["t"].to_numpy(dtype=float) * 10.0
        dt = np.maximum(minutes - row["onset_time"], 0.0)
        predicted = row["baseline_intensity"] + row["expression_amplitude"] * (
            np.exp(-row["protein_decay_rate"] * dt)
            - np.exp(-row["mrna_decay_rate"] * dt)
        )
        predicted = np.where(
            minutes < row["onset_time"], row["baseline_intensity"], predicted
        )
        ax.plot(minutes, predicted, color="green", alpha=0.1)
    ax.set_ylim(
        np.percentile(sc2["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc2["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc2_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[1, 0]
    for (pos, roi), trace in sc3.groupby(["pos", "roi"], sort=True):
        row = fit_df.loc[
            (fit_df["success"].astype(str).str.lower().eq("true"))
            & (fit_df["slide_channel"] == 3)
            & (fit_df["pos"] == pos)
            & (fit_df["roi"] == roi)
        ]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        minutes = trace["t"].to_numpy(dtype=float) * 10.0
        dt = np.maximum(minutes - row["onset_time"], 0.0)
        predicted = row["baseline_intensity"] + row["expression_amplitude"] * (
            np.exp(-row["protein_decay_rate"] * dt)
            - np.exp(-row["mrna_decay_rate"] * dt)
        )
        predicted = np.where(
            minutes < row["onset_time"], row["baseline_intensity"], predicted
        )
        ax.plot(minutes, predicted, color="green", alpha=0.1)
    ax.set_ylim(
        np.percentile(sc3["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc3["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc3_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[1, 1]
    for (pos, roi), trace in sc4.groupby(["pos", "roi"], sort=True):
        row = fit_df.loc[
            (fit_df["success"].astype(str).str.lower().eq("true"))
            & (fit_df["slide_channel"] == 4)
            & (fit_df["pos"] == pos)
            & (fit_df["roi"] == roi)
        ]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        minutes = trace["t"].to_numpy(dtype=float) * 10.0
        dt = np.maximum(minutes - row["onset_time"], 0.0)
        predicted = row["baseline_intensity"] + row["expression_amplitude"] * (
            np.exp(-row["protein_decay_rate"] * dt)
            - np.exp(-row["mrna_decay_rate"] * dt)
        )
        predicted = np.where(
            minutes < row["onset_time"], row["baseline_intensity"], predicted
        )
        ax.plot(minutes, predicted, color="green", alpha=0.1)
    ax.set_ylim(
        np.percentile(sc4["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc4["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc4_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    ax = axes[1, 2]
    for (pos, roi), trace in sc5.groupby(["pos", "roi"], sort=True):
        row = fit_df.loc[
            (fit_df["success"].astype(str).str.lower().eq("true"))
            & (fit_df["slide_channel"] == 5)
            & (fit_df["pos"] == pos)
            & (fit_df["roi"] == roi)
        ]
        if len(row) == 0:
            continue
        row = row.iloc[0]
        minutes = trace["t"].to_numpy(dtype=float) * 10.0
        dt = np.maximum(minutes - row["onset_time"], 0.0)
        predicted = row["baseline_intensity"] + row["expression_amplitude"] * (
            np.exp(-row["protein_decay_rate"] * dt)
            - np.exp(-row["mrna_decay_rate"] * dt)
        )
        predicted = np.where(
            minutes < row["onset_time"], row["baseline_intensity"], predicted
        )
        ax.plot(minutes, predicted, color="green", alpha=0.1)
    ax.set_ylim(
        np.percentile(sc5["corrected"].to_numpy(dtype=float), 5),
        np.percentile(sc5["corrected"].to_numpy(dtype=float), 95),
    )
    ax.set_title("sc5_ch1")
    ax.set_xlabel("minutes")
    ax.set_ylabel("corrected intensity")

    fig.tight_layout()
    fit_traces_fig = fig

    mo.vstack(
        [
            traces_fig,
            auc_linear_fig,
            auc_log_fig,
            protein_lifetime_fig,
            expression_rate_linear_fig,
            expression_rate_log_fig,
            fit_traces_fig,
        ]
    )


if __name__ == "__main__":
    app.run()
