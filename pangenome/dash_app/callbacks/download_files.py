import flask

from ..server import app


@app.server.route('/download_pangenome')
def download_json():
    return flask.send_file('../download/pangenome.json',
                           mimetype='text/csv',
                           attachment_filename='pangenome.json',
                           as_attachment=True)


@app.server.route('/download_csv')
def download_csv():
    return flask.send_file('../download/consensus.csv',
                           mimetype='text/csv',
                           attachment_filename='consensus.csv',
                           as_attachment=True)