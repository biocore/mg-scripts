import os


class HumanFilter:
    def __init__(self, seq_dir):
        self.seq_dir = seq_dir

    def run_param(self):
        for root, dirs, files in os.walk(self.seq_dir):
            for csv_file in files:
                # below is one method to search for a substring in a string
                if 'sav.csv' not in csv_file:
                    csv_file_path = os.path.join(root, csv_file)

    def process(self):
        pass
