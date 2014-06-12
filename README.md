

Installation
=======

### Dependencies

First (as this can run in parallel), download the datasets that the browser uses:

    wget http://broadinstitute.org/~bthomas/gencode.v19.annotation.gtf.gz
    wget http://broadinstitute.org/~bthomas/exac_chr20.vcf.gz

Follow these instructions to get Python and Homebrew installed on your Mac:
http://docs.python-guide.org/en/latest/starting/install/osx/

Install MongoDB:

    brew install mongodb

In a separate tab, start the mongo database server:

    mongod

This local server needs to be running at all times when you are working on the site.
You could do this in the background if you want or set up some startup service,
but I think it's easier just to open a tab you can monitor.

Finally, create a python virtual environment where the browser will live:

    virtualenv /path/to/env
    source /path/to/env/bin/activate

### Installation

First get the code:

    git clone https://github.com/brettpthomas/exac_browser.git
    cd exac_browser

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

For testing, you can open up an interactive shell with:

    python manage.py shell