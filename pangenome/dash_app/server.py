from flask import Flask
from dash import Dash

server = Flask('myproject')
app = Dash(server=server)