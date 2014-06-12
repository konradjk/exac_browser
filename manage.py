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


if __name__ == "__main__":
    manager.run()

