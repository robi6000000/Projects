import os
import sys
import pandas as pd
from sklearn.metrics import roc_auc_score, accuracy_score, confusion_matrix


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "Usage: python scripts/modelling/evaluate_single_predictions.py "
            "[prediction_csv] [output_csv]"
        )
        sys.exit(1)

    prediction_csv = sys.argv[1]
    output_csv = sys.argv[2]
    if not os.path.isfile(prediction_csv):
        print(f"Error: Prediction file '{prediction_csv}' does not exist.")
        sys.exit(1)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    df = pd.read_csv(prediction_csv)
    y_true = df['label'].values
    y_pred = df['probability'].values
    auc = roc_auc_score(y_true, y_pred)
    # save auc
    output_df = pd.DataFrame({
        'prediction_csv': [prediction_csv],
        'auc': [auc]
    })
    output_df.to_csv(output_csv, index=False)
    print(f"Evaluation complete for {prediction_csv}: AUC={auc:.4f}, results saved to {output_csv}")
