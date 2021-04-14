PROJECT = pysc

.PHONY: all test todo install develop doc version clean

all: version

install: version
	python setup.py install

user: version
	python setup.py install --user

test:
	python -m unittest discover -s tests

coverage:
	coverage run --source pysc -m unittest discover -s tests
	coverage report -m

clean:
	rm -rf ${PROJECT}/*pyc ${PROJECT}/*/*pyc ${PROJECT}/*/*so build
