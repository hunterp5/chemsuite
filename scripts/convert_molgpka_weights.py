#!/usr/bin/env python3
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

"""Convert MolGpKa's PyTorch weights to the numpy archives chemsuite ships.

MolGpKa (https://github.com/Xundrug/MolGpKa, MIT) saved its two models with torch's legacy
(pre-1.6, non-zip) format. This reads that format with the standard library and numpy, so the
conversion needs no torch install:

    python scripts/convert_molgpka_weights.py MolGpKa/models chemsuite/ionisation/data

It writes ``molgpka_acid.npz`` and ``molgpka_base.npz``. Only a developer runs this, once per
upstream weights change; the app reads the archives.
"""

from __future__ import annotations

import pickle
import struct
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

_DTYPES = {"FloatStorage": np.float32, "DoubleStorage": np.float64, "LongStorage": np.int64}


class _Storage:
    def __init__(self, kind: str) -> None:
        self.kind = kind
        self.data: np.ndarray | None = None


class _Tensor:
    def __init__(self, storage: _Storage, offset: int, size: tuple, stride: tuple) -> None:
        self.storage, self.offset, self.size, self.stride = storage, offset, size, stride

    def array(self) -> np.ndarray:
        data = self.storage.data
        assert data is not None
        if not self.size:
            return data[self.offset : self.offset + 1].reshape(()).copy()
        item = data.itemsize
        view = np.lib.stride_tricks.as_strided(
            data[self.offset :], shape=self.size, strides=[s * item for s in self.stride]
        )
        return np.ascontiguousarray(view)


def _rebuild_tensor(storage, offset, size, stride, *_rest):
    return _Tensor(storage, offset, tuple(size), tuple(stride))


class _Unpickler(pickle.Unpickler):
    def __init__(self, file, storages: dict[str, _Storage]) -> None:
        super().__init__(file)
        self.storages = storages

    def find_class(self, module: str, name: str):
        if module == "collections" and name == "OrderedDict":
            return OrderedDict
        if module == "torch._utils" and name in ("_rebuild_tensor_v2", "_rebuild_tensor"):
            return _rebuild_tensor
        if module == "torch" and name in _DTYPES:
            return name
        raise pickle.UnpicklingError(f"unexpected class {module}.{name} in the weights file")

    def persistent_load(self, saved_id):
        tag, kind, key, _location, _size, *_view = saved_id
        assert tag == "storage", saved_id
        return self.storages.setdefault(key, _Storage(kind))


def read_legacy_state_dict(path: Path) -> dict[str, np.ndarray]:
    """A ``torch.save``d state dict in the legacy format, as numpy arrays."""
    storages: dict[str, _Storage] = {}
    with path.open("rb") as file:
        # Every pickle goes through the restricted unpickler: it builds only plain data,
        # ordered dicts and tensors, never arbitrary classes.
        for _ in range(3):  # magic number, protocol version, system info
            _Unpickler(file, storages).load()
        state = _Unpickler(file, storages).load()
        keys = _Unpickler(file, storages).load()
        for key in keys:
            storage = storages[key]
            (count,) = struct.unpack("<q", file.read(8))
            dtype = np.dtype(_DTYPES[storage.kind]).newbyteorder("<")
            storage.data = np.frombuffer(file.read(count * dtype.itemsize), dtype=dtype)
    return {name: tensor.array() for name, tensor in state.items() if isinstance(tensor, _Tensor)}


def main(argv: list[str]) -> int:
    source, target = Path(argv[1]), Path(argv[2])
    target.mkdir(parents=True, exist_ok=True)
    for kind in ("acid", "base"):
        arrays = read_legacy_state_dict(source / f"weight_{kind}.pth")
        arrays = {k: v for k, v in arrays.items() if not k.endswith("num_batches_tracked")}
        out = target / f"molgpka_{kind}.npz"
        np.savez_compressed(out, **arrays)
        print(out, {k: v.shape for k, v in arrays.items()})
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
