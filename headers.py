#!/usr/bin/env python3
from __future__ import annotations

from Bio import SeqIO
import re

class header:
    def __init__(self, inp :str, out: str, pattern: list, new_id: str, new_desc: str) -> None:
        self.inp = inp
        self.out = out
        self.pattern = pattern
        self.new_id = new_id
        self. new_desc = new_desc

    @staticmethod
    def _renaming_logic(record: str, counter: int, switch_id: str, switch_description: str, outputfl: str) -> None:
        original_header = record.id
        original_description = record.description
        record.id = f'{switch_id}{counter}'
        record.description = f'{switch_description}{counter}'
        SeqIO.write(record, outputfl, 'fasta')
        print(f"\nHeader: {original_header} changed to {record.id}. Description changed from {original_description} to {record.description}.")

    def rename(self) -> None:
        """Temporary docstring: Uses regular expression to find pattern on headers."""

        with open(self.inp, "r") as fl, open(self.out, 'w') as outfl:
            headers = 0
            for line in fl:
                if line.startswith(">"):
                    headers += 1
            records1 = SeqIO.parse(self.inp, 'fasta')   # Use it to list id's.
            records2 = SeqIO.parse(self.inp, 'fasta')   # Use it for the rest of the operations.
            print(f"\nThere are {headers} headers in the input fasta.")
            iterate = 1
            record_ids = []
            for n in records1:
                x = n.id
                record_ids.append(x)

            invalid_patterns = []
            for record in records2:  # loop through records.
                for i in self.pattern:   # nested loop through specified patterns.
                    if not any(i in rec for rec in record_ids):
                        invalid_patterns.append(i)
                        continue

                    if iterate > 1: # only when the first pattern in the list is iterated.
                        if re.search(i, record.id):
                            self._renaming_logic(record = record, counter = iterate, switch_id = self.new_id, switch_description = self.new_desc, outputfl = outfl)
                            iterate += 1
                        else:
                            continue

                    else:   # first pattern iterable.
                        if re.search(i, record.id):
                            self._renaming_logic(record = record, counter = iterate, switch_id = self.new_id, switch_description = self.new_desc, outputfl = outfl)
                            iterate += 1
                        else:
                            continue

        # Final checks.
        if iterate == 1:
            print("\nNo headers with any of the defined patterns were found.\n")

        else:
            print(f"\n{iterate - 1} headers were modified and saved in {self.out}.\n")

        if len(invalid_patterns) > 0:
            str = ', '.join(set(invalid_patterns))
            print(f"The following specified patterns were not detected: {str}\n")


if __name__ == "__main__":
    inp_fl = "header_test_file.fasta"
    out_fl = "output.fasta"
    patterns = ["NC_045512.2", "header", "PPPPPPPPPP"]
    header(inp = inp_fl, out = out_fl, pattern = patterns, new_id = "foo", new_desc = "bar").rename()