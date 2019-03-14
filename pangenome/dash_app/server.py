from flask import Flask
from dash import Dash

server = Flask('pangenome')
app = Dash(server=server)