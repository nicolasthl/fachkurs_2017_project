"""

"""
import pandas
import sqlite3
import os.path

class ModelData:
    def __init__(self, datasource):
        if isinstance(datasource, str):
            if os.path.isdir(datasource):
                for file in os.listdir(datasource):
