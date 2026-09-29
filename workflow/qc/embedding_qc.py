"""Exploratory label prediction from already-fitted compound embeddings."""

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import LabelEncoder


def process_atc_code(code_string: str, level: int = 4) -> str:
    """Return the requested ATC hierarchy level (for example N03AX09 -> N03A)."""
    result = re.findall(r"[a-zA-Z]|[^a-zA-Z]+", code_string)
    return "".join(result[:level])


COUNT_THRESHOLDS = range(1, 50)
NUM_SPLITS = 5
SEED = 1234
np.random.seed(SEED)
ROOT = Path(__file__).resolve().parents[2]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--embeddings",
    type=Path,
    default=ROOT / "data/results/nextflow/molecule_embeddings.csv",
)
parser.add_argument(
    "--labels", type=Path, default=ROOT / "data/rawdata/HDD/colData.tsv"
)
parser.add_argument(
    "--outdir", type=Path, default=ROOT / "data/results/nextflow/qc"
)
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)

embs = pd.read_csv(args.embeddings, index_col=0)
if embs.empty or not embs.index.is_unique or embs.index.hasnans:
    raise ValueError(
        "Embeddings must contain unique, non-missing compound IDs and data"
    )
if not np.isfinite(embs.to_numpy(dtype=float)).all():
    raise ValueError("Embeddings contain non-finite values")
labels = pd.read_csv(
    args.labels,
    sep="\t",
    usecols=["HDD.Compound.ID", "Mechanism.of.Action", "ATC.Code"],
)
atc_labels = labels.dropna(subset=["ATC.Code"]).copy()
moa_labels = labels.dropna(subset=["Mechanism.of.Action"]).copy()
atc_labels["ATC.Code"] = atc_labels["ATC.Code"].map(
    lambda code: process_atc_code(code, level=3)
)
experiments = {
    "ATC": (atc_labels, "ATC.Code"),
    "MOA": (moa_labels, "Mechanism.of.Action"),
}
results = []
sizes = []
class_counts = []

for experiment, (data, label_col) in experiments.items():
    # Preserve the original cutoff definition: counts across the HDD metadata.
    counts = data[label_col].value_counts()
    for cutoff in COUNT_THRESHOLDS:
        eligible = data[data[label_col].isin(counts[counts >= cutoff].index)]
        eligible = eligible.set_index("HDD.Compound.ID")
        keep_mols = sorted(set(eligible.index).intersection(embs.index))
        eligible = eligible.loc[keep_mols, [label_col]]
        if not eligible.index.is_unique:
            raise ValueError(f"Duplicate compound labels for {experiment}")
        embedded_counts = eligible[label_col].value_counts()
        for label, count in embedded_counts.items():
            class_counts.append(
                {
                    "Target": experiment,
                    "Count Cutoff": cutoff,
                    "Label": label,
                    "HDD Count": int(counts[label]),
                    "Embedded Count": int(count),
                    "Present in Every Fold": bool(count >= NUM_SPLITS),
                }
            )
        # Preserve James's cohort: rare classes remain in the evaluation.
        # StratifiedKFold may warn when some classes have fewer than five samples.
        usable = eligible
        n_classes = usable[label_col].nunique()
        skip_reason = ""
        if n_classes < 2:
            skip_reason = "Fewer than two classes"
        elif len(usable) < NUM_SPLITS or embedded_counts.max() < NUM_SPLITS:
            skip_reason = "Insufficient samples for five-fold stratification"
        sizes.append(
            {
                "Target": experiment,
                "Num. Samples": len(usable),
                "Cutoff": cutoff,
                "Matched Samples": len(eligible),
                "Excluded Samples": len(eligible) - len(usable),
                "Num. Classes": n_classes,
                "CV Splits": NUM_SPLITS,
                "Status": "skipped" if skip_reason else "evaluated",
                "Reason": skip_reason,
            }
        )
        print(
            f"{experiment}, cutoff {cutoff}: {len(usable)} samples, {n_classes} classes"
            + (f"; skipped ({skip_reason})" if skip_reason else ""),
            flush=True,
        )
        if skip_reason:
            continue
        X = embs.loc[usable.index].to_numpy()
        y = LabelEncoder().fit_transform(usable[label_col])
        skf = StratifiedKFold(n_splits=NUM_SPLITS)
        models = {
            "RF": RandomForestClassifier(max_depth=2, random_state=0),
            "MLP": MLPClassifier(max_iter=1000),
            "KNN": KNeighborsClassifier(n_neighbors=10),
        }
        for fold, (train_index, test_index) in enumerate(
            skf.split(X, y), start=1
        ):
            for model_name, clf in models.items():
                clf.fit(X[train_index], y[train_index])
                preds = clf.predict(X[test_index])
                results.append(
                    {
                        "Target": experiment,
                        "Model": model_name,
                        "Count Cutoff": cutoff,
                        "Fold": fold,
                        "Bal. Acc": balanced_accuracy_score(
                            y[test_index], preds
                        ),
                        "Acc": accuracy_score(y[test_index], preds),
                    }
                )

results = pd.DataFrame(
    results,
    columns=["Target", "Model", "Count Cutoff", "Fold", "Bal. Acc", "Acc"],
)
sizes = pd.DataFrame(sizes)
results.to_csv(args.outdir / "classifier_performance.csv", index=False)
sizes.to_csv(args.outdir / "dataset_sizes.csv", index=False)
pd.DataFrame(
    class_counts,
    columns=[
        "Target",
        "Count Cutoff",
        "Label",
        "HDD Count",
        "Embedded Count",
        "Present in Every Fold",
    ],
).to_csv(args.outdir / "class_counts.csv", index=False)

for target in results["Target"].unique():
    target_results = results[results["Target"] == target]
    for metric, suffix, title in [
        ("Acc", "acc", "Accuracy"),
        ("Bal. Acc", "bal_acc", "Balanced accuracy"),
    ]:
        fig, ax = plt.subplots(figsize=(6.4, 4.8), layout="constrained")
        sns.lineplot(
            data=target_results,
            x="Count Cutoff",
            y=metric,
            hue="Model",
            errorbar=("ci", 95),
            seed=SEED,
            ax=ax,
        )
        ax.set(
            title=f"Cutoff vs. {title}\n{target} Prediction",
            xlabel="Count Cutoff",
        )
        fig.savefig(args.outdir / f"{target}_{suffix}.png", dpi=200)
        plt.close(fig)
print(f"Saved QC outputs to {args.outdir}")
