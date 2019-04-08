import json
from typing import Dict, Union

import jsonpickle
import pandas as pd
from pangenome.pang.fileformats.json.JSONPangenome import JSONPangenome


def unjsonify_jsonpangenome(jsonified_pangenome: str) -> JSONPangenome:
    return jsonpickle.decode(jsonified_pangenome)

def jsonify_dict(data: Dict) -> str:
    return json.dumps(data)

def unjsonify_dict(jsonified_data: str) -> Dict:
    return json.loads(jsonified_data)

def jsonify_df(df: pd.DataFrame) -> str:
    return df.to_json()

def unjsonify_df(jsonified_df: str) -> pd.DataFrame:
    return pd.read_json(jsonified_df)
