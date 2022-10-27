import argparse
from abc import abstractmethod, ABC
from dataclasses import dataclass
from typing import List, Optional, Dict


def _return_dual_topology_options(
    lines_a: List[str],
    lines_b: List[str],
    lines_ab: List[str],
    ifdef_a: str,
    ifdef_b: str,
    ifdef_ab: str,
):
    new_topology = [f"#ifdef {ifdef_a}\n"]
    new_topology.extend(lines_a)
    new_topology.append("#endif\n")
    new_topology.append(f"#ifdef {ifdef_b}\n")
    new_topology.extend(lines_b)
    new_topology.append("#endif\n")
    new_topology.append(f"#ifdef {ifdef_ab}\n")
    new_topology.extend(lines_ab)
    new_topology.append("#endif\n")

    return new_topology


def _is_comment_or_empty(line: str) -> bool:
    line = line.strip()
    return not line or line[0] == ";"


def _is_preprocessor_directive(line: str) -> bool:
    return line[0] == "#"


class LineBuffer:
    @abstractmethod
    def read_line(self, line: str) -> None:
        pass

    @abstractmethod
    def flush(self) -> List[str]:
        pass


class DefaultLineBuffer(LineBuffer):
    def __init__(self):
        self.buffer = []

    def read_line(self, line: str) -> None:
        self.buffer.append(line)

    def flush(self) -> List[str]:
        buffer = self.buffer
        self.buffer = []
        return buffer


class AtomTypeLineBuffer(LineBuffer):
    def __init__(self):
        self.buffer = [
            f"; {'name':<8} {'at.num':6} "
            f"{'mass':>10} {'charge':>10} "
            f"{'ptype':5} {'sigma':>11} {'epsilon':>11}\n"
        ]
        self.non_interacting = {}
        self.duplicates = {}
        self.extra_atoms = []
        self.extra_block = []

    def make_atom_optionally_non_interacting(
        self, name: str, alternative_name: str
    ) -> None:
        assert name not in self.non_interacting
        self.non_interacting[name] = alternative_name

    def duplicate_atom(self, name: str, alternative_name: str) -> None:
        assert name not in self.duplicates
        self.duplicates[name] = alternative_name

    def add_atom_type(
        self,
        name: str,
        atom_number: str,
        mass: float,
        charge: float,
        particle_type: str,
        sigma: float,
        epsilon: float,
    ):
        self.extra_atoms.append(
            f"{name:<10} {atom_number:>6} "
            f"{mass:10.6f} {charge:10.6f} "
            f"{particle_type:>5} {sigma:11.8f} {epsilon:11.8f}\n"
        )

    def add_follow_up_block(self, block: List[str]):
        self.extra_block.append("\n")
        self.extra_block.extend(block)

    def read_line(self, line: str) -> None:
        if _is_comment_or_empty(line):
            if len(self.buffer) > 1:
                # Skip empty lines or comments on top of block, otherwise keep them
                self.buffer.append(line)
            return
        if _is_preprocessor_directive(line):
            self.buffer.append(line)
            return

        line_list = line.split()
        if len(line_list) == 6:
            line_list.insert(1, "X")
        if len(line_list) == 8:
            line_list.pop(2)
        assert len(line_list) == 7

        (name, atom_number, mass, charge, particle_type, sigma, epsilon) = line.split()

        self.buffer.append(
            f"{name:<10} {atom_number:>6} "
            f"{float(mass):10.6f} {float(charge):10.6f} "
            f"{particle_type:>5} {float(sigma):11.8f} {float(epsilon):11.8f}\n"
        )

        if name in self.duplicates:
            self.buffer.append("; identical duplicate with different name\n")
            self.buffer.append(
                f"{self.duplicates[name]:<10} {atom_number:>6} "
                f"{float(mass):10.6f} {float(charge):10.6f} "
                f"{particle_type:>5} {float(sigma):11.8f} {float(epsilon):11.8f}\n"
            )
        if name in self.non_interacting:
            self.buffer.append("; non-interacting duplicate\n")
            self.buffer.append(
                f"{self.non_interacting[name]:<10} {atom_number:>6} "
                f"{float(mass):10.6f} {0.0:10.6f} "
                f"{particle_type:>5} {0.0:11.8f} {0.0:11.8f}\n"
            )

    def flush(self) -> List[str]:
        buffer = self.buffer
        while buffer[-1] == "\n":
            buffer.pop(-1)
        buffer.extend(self.extra_atoms)
        buffer.append("\n")
        self.buffer = []
        return buffer


class AtomLineBuffer(LineBuffer):
    def __init__(self):
        self.buffer_on = []
        self.buffer_off = []
        self.interaction_ifdef_name_on = None
        self.interaction_ifdef_name_off = None
        self.atom_type_switch = None

    def add_interaction_switch(self, ifdef_name_on: str, ifdef_name_off: str):
        self.interaction_ifdef_name_on = ifdef_name_on
        self.interaction_ifdef_name_off = ifdef_name_off

    def add_atom_type_switch(self, atom_type_switch: Dict[str, str]):
        self.atom_type_switch = atom_type_switch

    def read_line(self, line: str) -> None:
        if _is_comment_or_empty(line):
            return
        if _is_preprocessor_directive(line):
            self.buffer_on.append(line)
            self.buffer_off.append(line)
            return

        (
            nr,
            atype,
            resnr,
            residue,
            atom,
            cgnr,
            charge,
            mass,
        ) = line.split()

        atype_off = atype
        if self.atom_type_switch is not None and atype in self.atom_type_switch:
            atype_off = self.atom_type_switch[atype]

        self.buffer_on.append(
            f"{nr:>7} {atype:>14} "
            f"{resnr:>5} {residue:>7} {atom:>4} {cgnr:>5} "
            f"{float(charge):10.6f} {float(mass):7.4f}\n"
        )
        self.buffer_off.append(
            f"{nr:>7} {atype_off:>14} "
            f"{resnr:>5} {residue:>7} {atom:>4} {cgnr:>5} "
            f"{0.0:10.6f} {float(mass):7.4f}\n"
        )

    def flush(self) -> List[str]:
        buffer = [
            f";{'nr':>6} {'type':>14} "
            f"{'resnr':>5} {'residue':>7} {'atom':>4} {'cgnr':>5} "
            f"{'charge':>10} {'mass':>7}\n"
        ]

        closing_ifdef = []
        while _is_preprocessor_directive(self.buffer_on[-1]):
            closing_ifdef.insert(0, self.buffer_on.pop(-1))
            self.buffer_off.pop(-1)

        if self.interaction_ifdef_name_on is not None:
            buffer.append(f"#ifdef {self.interaction_ifdef_name_on}\n")
            buffer.extend(self.buffer_on)
            buffer.append("#endif\n")
            buffer.append(f"#ifdef {self.interaction_ifdef_name_off}\n")
            buffer.extend(self.buffer_off)
            buffer.append(f"#endif\n")
            buffer.append("\n")
            buffer.extend(closing_ifdef)
        else:
            buffer.extend(self.buffer_on)
        buffer.append("\n")
        self.buffer_on = []
        self.buffer_off = []
        return buffer


class DualLineBuffer(LineBuffer, ABC):
    def __init__(self, comment_header):
        self.buffer_a = []
        self.buffer_b = []
        self.buffer_ab = []
        self.comment_header = comment_header
        self._flush_line_by_line = True

    def _append_to_all(self, line):
        self.buffer_a.append(line)
        self.buffer_b.append(line)
        self.buffer_ab.append(line)

    def flush_line_by_line(self, value: bool) -> None:
        self._flush_line_by_line = value

    def flush(self) -> List[str]:
        length = len(self.buffer_a)
        assert (length == len(self.buffer_b)) and (length == len(self.buffer_ab))

        buffer = [self.comment_header]
        sub_buffer_a = []
        sub_buffer_b = []
        sub_buffer_ab = []

        if self._flush_line_by_line:
            for idx in range(length):
                if self.buffer_a[idx] == self.buffer_b[idx]:
                    # Line is identical in A and B
                    # Check if we non-empty sub-buffers, and flush them
                    if len(sub_buffer_a) > 0:
                        buffer.extend(
                            _return_dual_topology_options(
                                lines_a=sub_buffer_a,
                                lines_b=sub_buffer_b,
                                lines_ab=sub_buffer_ab,
                                ifdef_a="LIGAND_A",
                                ifdef_b="LIGAND_B",
                                ifdef_ab="LIGAND_AB",
                            )
                        )
                        sub_buffer_a.clear()
                        sub_buffer_b.clear()
                        sub_buffer_ab.clear()
                    # Add line directly to buffer
                    buffer.append(self.buffer_a[idx])
                else:
                    # Line is not identical, add line to sub-buffers
                    sub_buffer_a.append(self.buffer_a[idx])
                    sub_buffer_b.append(self.buffer_b[idx])
                    sub_buffer_ab.append(self.buffer_ab[idx])
        else:
            sub_buffer_a = self.buffer_a
            sub_buffer_b = self.buffer_b
            sub_buffer_ab = self.buffer_ab

        # Flush sub-buffers if needed
        if len(sub_buffer_a) > 0:
            buffer.extend(
                _return_dual_topology_options(
                    lines_a=sub_buffer_a,
                    lines_b=sub_buffer_b,
                    lines_ab=sub_buffer_ab,
                    ifdef_a="LIGAND_A",
                    ifdef_b="LIGAND_B",
                    ifdef_ab="LIGAND_AB",
                )
            )

        # Add empty line to finish block
        buffer.append("\n")

        # Empty buffers
        self.buffer_a = []
        self.buffer_b = []
        self.buffer_ab = []

        return buffer


class AtomDualLineBuffer(DualLineBuffer):
    def __init__(self):
        super().__init__(
            f";{'nr':>6} {'type':>14} "
            f"{'resnr':>5} {'residue':>7} {'atom':>4} {'cgnr':>5} "
            f"{'charge':>10} {'mass':>7} "
            f"{'typeB':>14} {'chargeB':>10} {'massB':>7}\n"
        )

    def read_line(self, line: str) -> None:
        if _is_comment_or_empty(line):
            return
        if _is_preprocessor_directive(line):
            self._append_to_all(line)
            return

        (
            nr,
            atype,
            resnr,
            residue,
            atom,
            cgnr,
            charge,
            mass,
            atype_b,
            charge_b,
            mass_b,
        ) = line.split()

        if (atype == atype_b) and (charge == charge_b) and (mass == mass_b):
            self._append_to_all(
                f"{nr:>7} {atype:>14} "
                f"{resnr:>5} {residue:>7} {atom:>4} {cgnr:>5} "
                f"{float(charge):10.6f} {float(mass):7.4f} "
                f"{atype_b:>14} {float(charge_b):10.6f} {float(mass_b):7.4f}\n"
            )
        else:
            self.buffer_a.append(
                f"{nr:>7} {atype:>14} "
                f"{resnr:>5} {residue:>7} {atom:>4} {cgnr:>5} "
                f"{float(charge):10.6f} {float(mass):7.4f} "
                f"{atype:>14} {float(charge):10.6f} {float(mass):7.4f}\n"
            )
            self.buffer_b.append(
                f"{nr:>7} {atype_b:>14} "
                f"{resnr:>5} {residue:>7} {atom:>4} {cgnr:>5} "
                f"{float(charge_b):10.6f} {float(mass_b):7.4f} "
                f"{atype_b:>14} {float(charge_b):10.6f} {float(mass_b):7.4f}\n"
            )
            self.buffer_ab.append(
                f"{nr:>7} {atype:>14} "
                f"{resnr:>5} {residue:>7} {atom:>4} {cgnr:>5} "
                f"{float(charge):10.6f} {float(mass):7.4f} "
                f"{atype_b:>14} {float(charge_b):10.6f} {float(mass_b):7.4f}\n"
            )


class BondDualLineBuffer(DualLineBuffer):
    def __init__(self):
        super().__init__(
            f";{'a1':>5} {'a2':>6} {'func':>6} "
            f"{'b_A':>14} {'k_A':>14} "
            f"{'b_B':>14} {'k_B':>14}\n"
        )

    def read_line(self, line: str) -> None:
        if _is_comment_or_empty(line):
            return
        if _is_preprocessor_directive(line):
            self._append_to_all(line)
            return

        (
            atom1,
            atom2,
            func,
            b_a,
            k_a,
            b_b,
            k_b,
        ) = line.split()
        if b_a != b_b or k_a != k_b:
            self.buffer_a.append(
                f"{atom1:>6} {atom2:>6} {func:>6} "
                f"{float(b_a):14.6f} {float(k_a):14.6f} "
                f"{float(b_a):14.6f} {float(k_a):14.6f}\n"
            )
            self.buffer_b.append(
                f"{atom1:>6} {atom2:>6} {func:>6} "
                f"{float(b_b):14.6f} {float(k_b):14.6f} "
                f"{float(b_b):14.6f} {float(k_b):14.6f}\n"
            )
            self.buffer_ab.append(
                f"{atom1:>6} {atom2:>6} {func:>6} "
                f"{float(b_a):14.6f} {float(k_a):14.6f} "
                f"{float(b_b):14.6f} {float(k_b):14.6f}\n"
            )
        else:
            self._append_to_all(
                f"{atom1:>6} {atom2:>6} {func:>6} "
                f"{float(b_a):14.6f} {float(k_a):14.6f} "
                f"{float(b_a):14.6f} {float(k_a):14.6f}\n"
            )


class AngleDualLineBuffer(DualLineBuffer):
    def __init__(self):
        super().__init__(
            f";{'a1':>5} {'a2':>6} {'a3':>6} {'func':>6} "
            f"{'theta_A':>14} {'k_A':>14} "
            f"{'theta_B':>14} {'k_B':>14}\n"
        )

    def read_line(self, line: str) -> None:
        if _is_comment_or_empty(line):
            return
        if _is_preprocessor_directive(line):
            self._append_to_all(line)
            return

        entries, comments = line.split(";")
        (
            atom1,
            atom2,
            atom3,
            func,
            theta_a,
            k_a,
            theta_b,
            k_b,
        ) = entries.split()
        if theta_a != theta_b or k_a != k_b:
            self.buffer_a.append(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {func:>6} "
                f"{float(theta_a):14.6f} {float(k_a):14.6f} "
                f"{float(theta_a):14.6f} {float(k_a):14.6f} "
                f"; {comments}"
            )
            self.buffer_b.append(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {func:>6} "
                f"{float(theta_b):14.6f} {float(k_b):14.6f} "
                f"{float(theta_b):14.6f} {float(k_b):14.6f} "
                f"; {comments}"
            )
            self.buffer_ab.append(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {func:>6} "
                f"{float(theta_a):14.6f} {float(k_a):14.6f} "
                f"{float(theta_b):14.6f} {float(k_b):14.6f} "
                f"; {comments}"
            )
        else:
            self._append_to_all(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {func:>6} "
                f"{float(theta_a):14.6f} {float(k_a):14.6f} "
                f"{float(theta_b):14.6f} {float(k_b):14.6f} "
                f";{comments}"
            )


class DihedralDualLineBuffer(DualLineBuffer):
    def __init__(self):
        super().__init__(
            f";{'a1':>5} {'a2':>6} {'a3':>6} {'a4':>6} {'func':>6} "
            f"{'phi_A':>5} {'k_A':>12} {'multi_A':>7} "
            f"{'phi_B':>5} {'k_B':>12} {'multi_B':>7}\n"
        )

    def read_line(self, line: str) -> None:
        if _is_comment_or_empty(line):
            return
        if _is_preprocessor_directive(line):
            self._append_to_all(line)
            return

        entries, comments = line.split(";")
        (
            atom1,
            atom2,
            atom3,
            atom4,
            func,
            phi_a,
            k_a,
            multi_a,
            phi_b,
            k_b,
            multi_b,
        ) = entries.split()

        if (phi_a == phi_b) and (k_a == k_b) and (multi_a == multi_b):
            self._append_to_all(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {atom4:>6} {func:>6} "
                f"{phi_a:>5} {float(k_a):12.7f} {multi_a:>7} "
                f"{phi_b:>5} {float(k_b):12.7f} {multi_b:>7} "
                f";{comments}"
            )
        else:
            self.buffer_a.append(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {atom4:>6} {func:>6} "
                f"{phi_a:>5} {float(k_a):12.7f} {multi_a:>7} "
                f"{phi_a:>5} {float(k_a):12.7f} {multi_a:>7} "
                f"; {comments}"
            )
            self.buffer_b.append(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {atom4:>6} {func:>6} "
                f"{phi_b:>5} {float(k_b):12.7f} {multi_b:>7} "
                f"{phi_b:>5} {float(k_b):12.7f} {multi_b:>7} "
                f"; {comments}"
            )
            self.buffer_ab.append(
                f"{atom1:>6} {atom2:>6} {atom3:>6} {atom4:>6} {func:>6} "
                f"{phi_a:>5} {float(k_a):12.7f} {multi_a:>7} "
                f"{phi_b:>5} {float(k_b):12.7f} {multi_b:>7} "
                f";{comments}"
            )


@dataclass
class BufferSwitch:
    previous_buffer_content: List[str]
    file_pointer_in_molecule_block: bool
    new_block: str
    new_buffer: LineBuffer


def _switch_buffer(
    current_block: Optional[str],
    new_block: str,
    buffer: LineBuffer,
    buffer_selection: Optional[Dict[str, LineBuffer]] = None,
    molecule_buffer_selection: Optional[Dict[str, LineBuffer]] = None,
    molecule_name: Optional[str] = None,
    file_pointer_is_in_molecule_block: bool = False,
) -> BufferSwitch:
    if molecule_buffer_selection is not None and molecule_name is None:
        raise ValueError(
            "Passed a molecule buffer selection, "
            "but no molecule name to match it against."
        )

    current_buffer = buffer.flush()
    if current_block == "moleculetype":
        relevant_lines = [
            line
            for line in current_buffer
            if (line.strip() and line.strip()[0] != ";" and line.strip()[0] != "[")
        ]
        assert len(relevant_lines) == 1
        molecule = relevant_lines[0].split()[0]
        file_pointer_is_in_molecule_block = molecule == molecule_name

    if new_block == "moleculetype":
        file_pointer_is_in_molecule_block = False

    if (
        file_pointer_is_in_molecule_block
        and molecule_buffer_selection is not None
        and new_block in molecule_buffer_selection
    ):
        new_buffer = molecule_buffer_selection[new_block]
    elif (
        not file_pointer_is_in_molecule_block
        and buffer_selection is not None
        and new_block in buffer_selection
    ):
        new_buffer = buffer_selection[new_block]
    else:
        new_buffer = DefaultLineBuffer()

    return BufferSwitch(
        previous_buffer_content=current_buffer,
        file_pointer_in_molecule_block=file_pointer_is_in_molecule_block,
        new_block=new_block,
        new_buffer=new_buffer,
    )


def _detect_block(line: str) -> Optional[str]:
    line = line.strip()
    if not line or line[0] != "[" or line[-1] != "]":
        return None
    return line[1:-2].strip()


def clean_dual_topology(
    input_topology: str, output_topology: str, ligand_name: str = "LIG"
) -> None:
    buffer = DefaultLineBuffer()
    topology = []
    current_block = None
    file_pointer_is_in_ligand_block = False

    buffer_selection = {
        "atoms": AtomDualLineBuffer(),
        "bonds": BondDualLineBuffer(),
        "angles": AngleDualLineBuffer(),
        "dihedrals": DihedralDualLineBuffer(),
    }
    buffer_selection["dihedrals"].flush_line_by_line(False)

    with open(input_topology) as in_file:
        for line in in_file:
            # Check if line is beginning of block
            block = _detect_block(line)

            if block is not None:
                # If yes, we found the beginning of a block and switch buffer
                buffer_switch = _switch_buffer(
                    current_block=current_block,
                    new_block=block,
                    buffer=buffer,
                    molecule_buffer_selection=buffer_selection,
                    molecule_name=ligand_name,
                    file_pointer_is_in_molecule_block=file_pointer_is_in_ligand_block,
                )
                topology.extend(buffer_switch.previous_buffer_content)
                topology.append(line)
                file_pointer_is_in_ligand_block = (
                    buffer_switch.file_pointer_in_molecule_block
                )
                current_block = buffer_switch.new_block
                buffer = buffer_switch.new_buffer
            else:
                # Otherwise, simply read line with current buffer
                buffer.read_line(line)

    # Flush last buffer
    topology.extend(buffer.flush())

    # Write new topology file
    with open(output_topology, "w") as out_file:
        for line in topology:
            out_file.write(line)


def duplicate_water_molecule_and_add_atom_types(
    input_topology: str,
    output_topology: str,
    water_name: str = "SOL",
    special_water_name: str = "MOL",
    atom_duplication: Optional[Dict[str, str]] = None,
    atom_optionally_non_interacting: Optional[Dict[str, str]] = None,
    water_atom_type_changes: Optional[Dict[str, str]] = None,
) -> None:
    buffer = DefaultLineBuffer()
    topology = []
    current_block = None
    file_pointer_is_in_water_block = False
    buffer_selection = {"atomtypes": AtomTypeLineBuffer()}
    for atom_type in atom_duplication:
        buffer_selection["atomtypes"].duplicate_atom(
            atom_type, atom_duplication[atom_type]
        )
    for atom_type in atom_optionally_non_interacting:
        buffer_selection["atomtypes"].make_atom_optionally_non_interacting(
            atom_type, atom_optionally_non_interacting[atom_type]
        )
    buffer_selection["atomtypes"].add_atom_type(
        name="attach",
        atom_number="AT",
        mass=0.0,
        charge=0.0,
        particle_type="D",
        sigma=0.0,
        epsilon=0.0,
    )
    buffer_selection["atomtypes"].add_atom_type(
        name="attachX",
        atom_number="AT",
        mass=0.0,
        charge=0.0,
        particle_type="D",
        sigma=0.0,
        epsilon=0.0,
    )
    buffer_selection["atomtypes"].add_follow_up_block(
        [
            "[ nonbond_params ]\n",
            "; i      j         func       V(sigma)     W(eps)\n",
            "attachX  OWX       1          -0.315365    0.648520\n",
        ]
    )

    with open(input_topology) as in_file:
        water_topology = []
        for line in in_file:
            # Check if line is beginning of block
            block = _detect_block(line)

            if block is not None:
                # If yes, we found the beginning of a block and switch buffer
                buffer_switch = _switch_buffer(
                    current_block=current_block,
                    new_block=block,
                    buffer=buffer,
                    buffer_selection=buffer_selection,
                    molecule_name=water_name,
                    file_pointer_is_in_molecule_block=file_pointer_is_in_water_block,
                )
                if (
                    not file_pointer_is_in_water_block
                    and buffer_switch.file_pointer_in_molecule_block
                ):
                    # This is the header block of the water, we need to retrieve its
                    # tag from the default topology buffer
                    water_topology.append(topology.pop(-1))
                if (
                    file_pointer_is_in_water_block
                    or buffer_switch.file_pointer_in_molecule_block
                ):
                    # We were in water block before, so the previous buffer is part of
                    # it, or we are just starting it, which means that the current
                    # buffer contains the header
                    water_topology.extend(buffer_switch.previous_buffer_content)
                else:
                    # This is not part of the water block, we can just add it to the
                    # regular topology
                    topology.extend(buffer_switch.previous_buffer_content)

                # Check if we are in the water block now
                file_pointer_is_in_water_block = (
                    buffer_switch.file_pointer_in_molecule_block
                )
                if file_pointer_is_in_water_block:
                    # We are in water block now, add line to water topology
                    water_topology.append(line)
                else:
                    # We are not in water block
                    if water_topology:
                        # We just finished the water block, we need to transfer its
                        # edited version to the regular topology
                        # First write new special water block
                        for water_line in water_topology:
                            if water_name in water_line:
                                water_line = water_line.replace(
                                    water_name, special_water_name
                                )
                            for atom_type in water_atom_type_changes:
                                if atom_type in water_line:
                                    water_line = water_line.replace(
                                        atom_type, water_atom_type_changes[atom_type]
                                    )
                            topology.append(water_line)
                        # Then write normal water block
                        topology.extend(water_topology)
                        # Empty water buffer
                        water_topology = []
                    topology.append(line)

                current_block = buffer_switch.new_block
                buffer = buffer_switch.new_buffer

                if block == "atomtypes" and "atomtypes" in buffer_selection:
                    # This block might be there repeatedly, we only want to add things
                    # to the first. This is assuming that the atomtypes block describing
                    # the protein and the water comes first. Ligand atom types need to
                    # be either described in the same block or in a later block.
                    # TODO: If things are in a different order, this doesn't work!
                    buffer_selection.pop("atomtypes")
            else:
                # Otherwise, simply read line with current buffer
                buffer.read_line(line)

    # Flush last buffer
    topology.extend(buffer.flush())

    # Write new topology file
    with open(output_topology, "w") as out_file:
        for line in topology:
            out_file.write(line)


def add_non_interacting_water(
    input_topology: str,
    output_topology: str,
    special_water_name: str = "MOL",
    water_atom_type_changes: Optional[Dict[str, str]] = None,
) -> None:
    buffer = DefaultLineBuffer()
    topology = []
    current_block = None
    file_pointer_is_in_water_block = False
    molecule_buffer_selection = {"atoms": AtomLineBuffer()}
    molecule_buffer_selection["atoms"].add_interaction_switch(
        ifdef_name_on=f"{special_water_name}_INTERACTING",
        ifdef_name_off=f"{special_water_name}_NONINTERACTING",
    )
    molecule_buffer_selection["atoms"].add_atom_type_switch(water_atom_type_changes)

    with open(input_topology) as in_file:
        for line in in_file:
            # Check if line is beginning of block
            block = _detect_block(line)
            if block is not None:
                # If yes, we found the beginning of a block and switch buffer
                buffer_switch = _switch_buffer(
                    current_block=current_block,
                    new_block=block,
                    buffer=buffer,
                    molecule_buffer_selection=molecule_buffer_selection,
                    molecule_name=special_water_name,
                    file_pointer_is_in_molecule_block=file_pointer_is_in_water_block,
                )
                topology.extend(buffer_switch.previous_buffer_content)
                topology.append(line)
                file_pointer_is_in_water_block = (
                    buffer_switch.file_pointer_in_molecule_block
                )
                current_block = buffer_switch.new_block
                buffer = buffer_switch.new_buffer
            else:
                # Otherwise, simply read line with current buffer
                buffer.read_line(line)

    # Flush last buffer
    topology.extend(buffer.flush())

    # Write new topology file
    with open(output_topology, "w") as out_file:
        for line in topology:
            out_file.write(line)


def command_line_entry_point():
    parser = argparse.ArgumentParser(
        description="The water NES workflow needs fine-grained control over the "
        "lambda parameters (e.g., we want to use bonded interactions to "
        "restrain the water molecule, and need to be able to modify that "
        "specific interaction independently of other lambda alchemical "
        "transformations). "
        "It is therefore necessary to have three definitions of the dual "
        "topology: One with only state A, one with only state B, and one "
        "with both states, all guarded with #ifdef statements. "
        "This script takes a topology with a ligand having A and B "
        "parameters (as, for example, prepared by pmx) as input and saves "
        "a new file with the flexibility required by the water NES workflow."
    )
    parser.add_argument(
        "input_topology_file", type=str, help="The file to read the topology from"
    )
    parser.add_argument(
        "output_topology_file", type=str, help="The file to write the topology to"
    )
    parser.add_argument(
        "--ligand_name",
        type=str,
        metavar="LIG",
        default="LIG",
        help="The name of the ligand molecule in the input topology. Default: LIG",
    )
    args = parser.parse_args()
    clean_dual_topology(
        input_topology=args.input_topology_file,
        output_topology=args.output_topology_file + ".1",
        ligand_name=args.ligand_name,
    )
    duplicate_water_molecule_and_add_atom_types(
        input_topology=args.output_topology_file + ".1",
        output_topology=args.output_topology_file + ".2",
        water_name="SOL",
        special_water_name="MOL",
        atom_duplication={"OW": "OWX"},
        atom_optionally_non_interacting={"OW": "OWD"},
        water_atom_type_changes={"OW": "OWX"},
    )
    add_non_interacting_water(
        input_topology=args.output_topology_file + ".2",
        output_topology=args.output_topology_file,
        special_water_name="MOL",
        water_atom_type_changes={"OWX": "OWD"},
    )


if __name__ == "__main__":
    command_line_entry_point()
