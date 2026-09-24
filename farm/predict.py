"""Top-level convenience entry point: ``from farm import predict``."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .model import load_model


def predict(
    data,
    model: str = "combined",
    threshold: float | None = None,
    return_frame: bool = True,
):
    """Score variants for resistance association.

    Parameters
    ----------
    data
        A :class:`pandas.DataFrame` carrying the 25 FARM feature columns, or a path to
        a CSV holding them. Identifier columns such as ``gene`` and
        ``mutation_oneletter`` are carried through to the output when present.
    model
        ``"combined"`` (default, the deployment model), ``"essential"`` or
        ``"nonessential"``.
    threshold
        Override the model's calibrated decision threshold. Leave as ``None`` to use
        the manuscript value.
    return_frame
        When ``True`` return a frame with scores and calls; when ``False`` return the
        raw probability array.

    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        With ``return_frame=True``: columns ``resistance_score`` (probability),
        ``predicted_resistant`` (0/1), ``threshold`` and ``farm_model``, plus any
        identifier columns found on the input.

    Examples
    --------
    >>> from farm import predict
    >>> scored = predict("variants_with_features.csv")           # doctest: +SKIP
    >>> scored.sort_values("resistance_score", ascending=False)  # doctest: +SKIP
    """
    if isinstance(data, (str, Path)):
        data = pd.read_csv(data)

    fm = load_model(model)
    if return_frame:
        return fm.predict_frame(data, threshold=threshold)
    return fm.predict_proba(data)
