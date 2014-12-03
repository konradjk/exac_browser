from flask.ext.script import Manager
from exac import app
import exac

manager = Manager(app)


@manager.command
def hello():
    print "hello"


@manager.command
def load_db():
    exac.load_db()


@manager.command
def load_base_coverage():
    exac.load_base_coverage()


@manager.command
def load_variants_file():
    exac.load_variants_file()


@manager.command
def load_gene_models():
    exac.load_gene_models()


@manager.command
def load_dbsnp():
    exac.load_dbsnp()


@manager.command
def create_cache():
    exac.create_cache()


if __name__ == "__main__":
    manager.run()

