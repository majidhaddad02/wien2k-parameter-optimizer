"""
WIEN2k case.struct Parser with NN computation.

References:
  - WIEN2k User's Guide, struct file format
"""

import math
import re
from dataclasses import dataclass, field


@dataclass
class Atom:
    index: int
    element: str
    z: int
    mult: int
    rmt: float
    position: list
    isplit: int = 0


@dataclass
class Structure:
    title: str = ""
    lattice_type: str = ""
    mode: str = ""
    a: float = 0.0
    b: float = 0.0
    c: float = 0.0
    alpha: float = 90.0
    beta: float = 90.0
    gamma: float = 90.0
    atoms: list = field(default_factory=list)
    spacegroup: int = 1

    @property
    def lattice_vectors(self):
        al = math.radians(self.alpha)
        be = math.radians(self.beta)
        ga = math.radians(self.gamma)
        ca, cb, cg = math.cos(al), math.cos(be), math.cos(ga)
        sg = math.sin(ga)

        a_vec = [self.a, 0.0, 0.0]
        b_vec = [self.b * cg, self.b * sg, 0.0]

        if sg == 0:
            c_vec = [0.0, 0.0, self.c]
        else:
            c1 = self.c * cb
            c2 = self.c * (ca - cb * cg) / sg
            v = 1.0 - ca**2 - cb**2 - cg**2 + 2 * ca * cb * cg
            c3 = self.c * math.sqrt(max(v, 0.0)) / sg
            c_vec = [c1, c2, c3]
        return a_vec, b_vec, c_vec

    @property
    def volume(self):
        a_vec, b_vec, c_vec = self.lattice_vectors
        return abs(
            a_vec[0] * (b_vec[1] * c_vec[2] - b_vec[2] * c_vec[1])
            - a_vec[1] * (b_vec[0] * c_vec[2] - b_vec[2] * c_vec[0])
            + a_vec[2] * (b_vec[0] * c_vec[1] - b_vec[1] * c_vec[0])
        )

    @property
    def reciprocal_vectors(self):
        a_vec, b_vec, c_vec = self.lattice_vectors
        V = self.volume
        if V == 0:
            return [1, 0, 0], [0, 1, 0], [0, 0, 1]
        b1 = [
            2 * math.pi * (b_vec[1] * c_vec[2] - b_vec[2] * c_vec[1]) / V,
            2 * math.pi * (b_vec[2] * c_vec[0] - b_vec[0] * c_vec[2]) / V,
            2 * math.pi * (b_vec[0] * c_vec[1] - b_vec[1] * c_vec[0]) / V,
        ]
        b2 = [
            2 * math.pi * (c_vec[1] * a_vec[2] - c_vec[2] * a_vec[1]) / V,
            2 * math.pi * (c_vec[2] * a_vec[0] - c_vec[0] * a_vec[2]) / V,
            2 * math.pi * (c_vec[0] * a_vec[1] - c_vec[1] * a_vec[0]) / V,
        ]
        b3 = [
            2 * math.pi * (a_vec[1] * b_vec[2] - a_vec[2] * b_vec[1]) / V,
            2 * math.pi * (a_vec[2] * b_vec[0] - a_vec[0] * b_vec[2]) / V,
            2 * math.pi * (a_vec[0] * b_vec[1] - a_vec[1] * b_vec[0]) / V,
        ]
        return b1, b2, b3

    @property
    def reciprocal_lengths(self):
        b1, b2, b3 = self.reciprocal_vectors
        def _norm(v):
            return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        return _norm(b1), _norm(b2), _norm(b3)

    @property
    def bz_volume(self):
        return (2 * math.pi)**3 / max(self.volume, 1e-10)

    @property
    def num_atoms_primitive(self):
        return sum(a.mult for a in self.atoms)


def parse_struct(filepath: str) -> Structure:
    with open(filepath, "r") as f:
        lines = [l.rstrip("\n").rstrip("\r") for l in f]

    s = Structure()
    idx = 0
    s.title = lines[idx].strip()
    idx += 1

    s.lattice_type = lines[idx].strip()
    idx += 1

    while idx < len(lines) and lines[idx].strip() == "":
        idx += 1
    if idx < len(lines):
        s.mode = lines[idx].strip()
        idx += 1

    while idx < len(lines) and lines[idx].strip() == "":
        idx += 1
    if idx < len(lines):
        parts = lines[idx].split()
        nums = []
        for p in parts:
            try:
                nums.append(float(p))
            except ValueError:
                nums.append(0.0)
        if len(nums) >= 6:
            s.a, s.b, s.c = nums[0], nums[1], nums[2]
            s.alpha, s.beta, s.gamma = nums[3], nums[4], nums[5]
        idx += 1

    while idx < len(lines):
        line = lines[idx].strip()
        idx += 1
        if not line or not line.startswith("ATOM"):
            continue

        pos_match = re.search(
            r"X\s*=\s*([0-9.eE+\-]+)\s+Y\s*=\s*([0-9.eE+\-]+)\s+Z\s*=\s*([0-9.eE+\-]+)",
            line
        )
        pos = [0.0, 0.0, 0.0]
        if pos_match:
            pos = [float(pos_match.group(i)) for i in (1, 2, 3)]

        mult, isplit = 1, 0
        if idx < len(lines):
            ml = lines[idx].strip()
            mm = re.search(r"MULT\s*=\s*(\d+)", ml)
            im = re.search(r"ISPLIT\s*=\s*(\d+)", ml)
            if mm:
                mult = int(mm.group(1))
            if im:
                isplit = int(im.group(1))
            idx += 1

        element = ""
        z_val = 0
        rmt_val = 2.0
        if idx < len(lines):
            il = lines[idx].strip()
            ip = il.split()
            if ip:
                raw = re.sub(r'[\d_]', '', ip[0])
                element = raw.strip().title()
            for ti, tok in enumerate(ip):
                tu = tok.upper()
                if tu == "Z:" and ti + 1 < len(ip):
                    try:
                        z_val = int(float(ip[ti + 1].rstrip(",")))
                    except (ValueError, IndexError):
                        pass
                elif tu.startswith("Z:") and tu != "Z:":
                    try:
                        z_val = int(float(tu.split(":")[1]))
                    except (ValueError, IndexError):
                        pass
                if tu == "RMT=" and ti + 1 < len(ip):
                    try:
                        rmt_val = float(ip[ti + 1])
                    except (ValueError, IndexError):
                        pass
                elif tu.startswith("RMT=") and tu != "RMT=":
                    try:
                        rmt_val = float(tu.split("=")[1])
                    except (ValueError, IndexError):
                        pass
            idx += 1

        while idx < len(lines):
            la = lines[idx].strip()
            if la.startswith("ATOM") or not la:
                break
            idx += 1

        if not element:
            continue

        atom = Atom(
            index=len(s.atoms),
            element=element,
            z=z_val,
            mult=mult,
            rmt=rmt_val,
            position=pos,
            isplit=isplit,
        )
        s.atoms.append(atom)

    return s


def compute_pairwise_min_distances(structure: Structure):
    """Min distance between every pair of atoms including [-2,-1,0,1,2] cell images."""
    import math as _math
    a_vec, b_vec, c_vec = structure.lattice_vectors
    n = len(structure.atoms)
    pairwise = {}

    for i in range(n):
        ri = structure.atoms[i].position
        for j in range(i + 1, n):
            rj = structure.atoms[j].position
            md = float("inf")
            for na in (-2, -1, 0, 1, 2):
                for nb in (-2, -1, 0, 1, 2):
                    for nc in (-2, -1, 0, 1, 2):
                        dx = (ri[0] - rj[0] - na) * a_vec[0] \
                           + (ri[1] - rj[1] - nb) * b_vec[0] \
                           + (ri[2] - rj[2] - nc) * c_vec[0]
                        dy = (ri[0] - rj[0] - na) * a_vec[1] \
                           + (ri[1] - rj[1] - nb) * b_vec[1] \
                           + (ri[2] - rj[2] - nc) * c_vec[1]
                        dz = (ri[0] - rj[0] - na) * a_vec[2] \
                           + (ri[1] - rj[1] - nb) * b_vec[2] \
                           + (ri[2] - rj[2] - nc) * c_vec[2]
                        d = _math.sqrt(dx*dx + dy*dy + dz*dz)
                        if d < md:
                            md = d
            pairwise[(i, j)] = md
            pairwise[(j, i)] = md
    return pairwise
