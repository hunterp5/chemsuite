# This file is part of chemsuite.
# Copyright (C) 2026 Hunter Picard
#
# chemsuite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# chemsuite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with chemsuite.  If not, see <https://www.gnu.org/licenses/>.

"""The Uni-pKa network's inference in numpy: a microstate's free energy from its 3D atoms.

A port of the forward pass of Uni-pKa's ``unimol_pka`` model (Luo et al., JACS Au 2024,
doi:10.1021/jacsau.4c00271; https://github.com/dptech-corp/Uni-pKa, Apache-2.0), which is
Uni-Mol's encoder (Zhou et al., ICLR 2023; https://github.com/deepmodeling/Uni-Mol, MIT) with an
embedding of formal charges. chemsuite trains its own weights (docs/UNIPKA_PLAN.md); they are
converted from PyTorch with ``scripts/unipka/convert_weights.py``, so the app needs no torch.

The input is one microstate, hydrogens explicit, in one or more conformers:

* tokens: each atom's element in :data:`ATOMS` (``[UNK]`` for others), and its formal charge
  in :data:`CHARGES`, framed by ``[CLS]`` and ``[SEP]``;
* coordinates centred on the atoms' mean, ``[CLS]`` and ``[SEP]`` at the origin;
* the pair bias: each distance through a Gaussian basis whose scale and shift depend on the
  pair's two element tokens, projected to one bias per attention head.

Each of the pre-LN layers adds its attention logits (the pair bias included) to the next
layer's bias. The ``[CLS]`` row gives the free energy in kT; conformers are averaged.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.special import erf

#: Uni-Mol's dictionary (``dict.txt``) and Uni-pKa's charge dictionary (``dict_charge.txt``)
ATOMS = ("[PAD]", "[CLS]", "[SEP]", "[UNK]", "C", "N", "O", "S", "H", "Cl", "F", "Br", "I", "Si",
         "P", "B", "Na", "K", "Al", "Ca", "Sn", "As", "Hg", "Fe", "Zn", "Cr", "Se", "Gd", "Au",
         "Li", "[MASK]")  # fmt: skip
CHARGES = ("[PAD]", "[CLS]", "[SEP]", "[UNK]", "1", "0", "-1", "[MASK]")
_CLS, _SEP, _UNK = 1, 2, 3
_EPS = 1e-5  # Uni-Core's LayerNorm
_SQRT_2PI = (2 * 3.14159) ** 0.5  # upstream's constant, π to five figures


@dataclass(frozen=True)
class Network:
    """One fold's parameters, by their PyTorch names less the ``unimol.`` prefix."""

    p: dict[str, np.ndarray]
    layers: int
    heads: int
    head: str  # the regression head's name

    @classmethod
    def load(cls, path: Path, fold: int | None = None) -> Network:
        """The fold *fold* of the archive at *path* (keys ``f<k>/<name>``), or its only one."""
        with np.load(path) as archive:
            folds = sorted({k.split("/", 1)[0] for k in archive.files if k.startswith("f")})
            prefix = folds[0] if fold is None else f"f{fold}"
            p = {k.split("/", 1)[1]: archive[k].astype(np.float32)
                 for k in archive.files if k.startswith(prefix + "/")}  # fmt: skip
        layers = 1 + max(int(k.split(".")[2]) for k in p if k.startswith("encoder.layers."))
        heads = p["gbf_proj.linear2.weight"].shape[0]
        (head,) = {k.split(".")[1] for k in p if k.startswith("classification_heads.")}
        return cls(p, layers, heads, head)


def tokens(atoms: list[str], charges: list[int]) -> tuple[np.ndarray, np.ndarray]:
    """The element and charge tokens of a microstate, ``[CLS]`` first and ``[SEP]`` last."""
    index = {a: k for k, a in enumerate(ATOMS)}
    charge_index = {c: k for k, c in enumerate(CHARGES)}
    atom_tokens = [_CLS] + [index.get(a, _UNK) for a in atoms] + [_SEP]
    charge_tokens = [_CLS] + [charge_index.get(str(q), _UNK) for q in charges] + [_SEP]
    return np.array(atom_tokens), np.array(charge_tokens)


def _layer_norm(x: np.ndarray, w: np.ndarray, b: np.ndarray) -> np.ndarray:
    mean = x.mean(-1, keepdims=True)
    var = ((x - mean) ** 2).mean(-1, keepdims=True)
    return (x - mean) / np.sqrt(var + _EPS) * w + b


def _gelu(x: np.ndarray) -> np.ndarray:
    return 0.5 * x * (1.0 + erf(x / np.sqrt(2.0)))


def _linear(x: np.ndarray, p: dict[str, np.ndarray], name: str) -> np.ndarray:
    return x @ p[f"{name}.weight"].T + p[f"{name}.bias"]


def _softmax(x: np.ndarray) -> np.ndarray:
    x = x - x.max(-1, keepdims=True)
    e = np.exp(x)
    return e / e.sum(-1, keepdims=True)


def free_energy(
    net: Network, atoms: list[str], charges: list[int], coords: np.ndarray
) -> np.ndarray:
    """The free energy (kT) of a microstate in each of its conformers *coords* (C, n, 3)."""
    p = net.p
    tok, chg = tokens(atoms, charges)
    c = coords.shape[0]
    xyz = coords - coords.mean(axis=1, keepdims=True)
    xyz = np.concatenate([np.zeros((c, 1, 3)), xyz, np.zeros((c, 1, 3))], axis=1)
    n = xyz.shape[1]
    dist = np.linalg.norm(xyz[:, :, None] - xyz[:, None], axis=-1).astype(np.float32)
    # the pair bias: Gaussian basis on each distance, scaled and shifted by the pair's elements
    edge = tok[:, None] * len(ATOMS) + tok[None, :]
    mul = p["gbf.mul.weight"][edge, 0]
    shift = p["gbf.bias.weight"][edge, 0]
    mean = p["gbf.means.weight"].reshape(-1)
    std = np.abs(p["gbf.stds.weight"].reshape(-1)) + 1e-5
    scaled = (mul * dist + shift)[..., None]
    basis = np.exp(-0.5 * ((scaled - mean) / std) ** 2) / (_SQRT_2PI * std)
    bias = _linear(_gelu(_linear(basis, p, "gbf_proj.linear1")), p, "gbf_proj.linear2")
    bias = bias.transpose(0, 3, 1, 2)  # (C, heads, n, n)
    # the encoder
    x = p["embed_tokens.weight"][tok] + p["embed_charges.weight"][chg]
    x = np.repeat(x[None], c, axis=0)  # contiguous: a broadcast view makes numpy skip BLAS
    x = _layer_norm(x, p["encoder.emb_layer_norm.weight"], p["encoder.emb_layer_norm.bias"])
    h = net.heads
    d = x.shape[-1] // h
    for k in range(net.layers):
        name = f"encoder.layers.{k}"
        y = _layer_norm(x, p[f"{name}.self_attn_layer_norm.weight"],
                        p[f"{name}.self_attn_layer_norm.bias"])  # fmt: skip
        q, kk, v = np.split(_linear(y, p, f"{name}.self_attn.in_proj"), 3, axis=-1)
        q = q.reshape(c, n, h, d).transpose(0, 2, 1, 3) * d**-0.5
        kk = kk.reshape(c, n, h, d).transpose(0, 2, 1, 3)
        v = v.reshape(c, n, h, d).transpose(0, 2, 1, 3)
        bias = q @ kk.transpose(0, 1, 3, 2) + bias  # the logits become the next layer's bias
        o = (_softmax(bias) @ v).transpose(0, 2, 1, 3).reshape(c, n, h * d)
        x = x + _linear(o, p, f"{name}.self_attn.out_proj")
        y = _layer_norm(x, p[f"{name}.final_layer_norm.weight"], p[f"{name}.final_layer_norm.bias"])
        x = x + _linear(_gelu(_linear(y, p, f"{name}.fc1")), p, f"{name}.fc2")
    x = _layer_norm(x, p["encoder.final_layer_norm.weight"], p["encoder.final_layer_norm.bias"])
    head = f"classification_heads.{net.head}"
    out = _linear(np.tanh(_linear(x[:, 0], p, f"{head}.dense")), p, f"{head}.out_proj")
    return out[:, 0].astype(np.float64)
