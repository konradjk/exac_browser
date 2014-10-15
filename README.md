

Installation
=======

### Getting the code

Create a directory to put all this stuff in. This will serve as the parent directory of the actual exac_browser repository 

    mkdir exac
    cd exac

First (as this can run in parallel), get the datasets that the browser uses:

    scp tin:/humgen/atgu1/fs03/konradk/exac_browser/exac_browser.tar.gz .
    tar zxvf exac_browser.tar.gz

Now clone the repo: 

    git clone https://github.com/brettpthomas/exac_browser.git

### Dependencies

Follow these instructions to get Python and Homebrew installed on your Mac:
http://docs.python-guide.org/en/latest/starting/install/osx/

Install MongoDB:

    brew install mongodb
    # or
    sudo port install mongodb

Create a directory to hold your mongo database files: 

    mkdir database

In a separate tab, start the mongo database server:

    mongod --dbpath database

This local server needs to be running at all times when you are working on the site.
You could do this in the background if you want or set up some startup service,
but I think it's easier just to open a tab you can monitor.

Finally, you may want to keep the system in a virtualenv:

    sudo port install py27-virtualenv # Or whatever version

If so, you can create a python virtual environment where the browser will live:

    mkdir exac_env
    virtualenv exac_env
    source exac_env/bin/activate

### Installation

Install the python requirements:

    pip install -r requirements.txt

Note that this installs xBrowse too.

### Setup

First, make sure those data files you downloaded earlier are in the same directory as `exac.py`.

At this point, it's also probably worth quickly checking out the code structure if you haven't already :)

Now we must load the database from those flat files.
This is a single command.
It can take a while, but there are some cute progress bars to keep you company:

    python manage.py load_db
    # TODO: ./manage.py doesn't work for some reason - I guess numbered args are used somewhere.

You won't have to run this often - most changes won't require rebuilding the database.
That said, this is (and will remain) idempotent,
so you can run it again at any time if you think something might be wrong - it will reload the database from scratch.

### Running the site

Note that if you are revisiting the site after a break, make sure your virtualenv is `activate`'d.

You can run the development server with:

    python exac.py

And visit on your browser:

    http://localhost:5000
    http://localhost:5000/gene/ENSG00000237683
    http://localhost:5000/variant/20-76735-A-T


For testing, you can open up an interactive shell with:

    python manage.py shell