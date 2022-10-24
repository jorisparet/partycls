PROJECT = partycls

.PHONY: all fortran test install coverage version clean

all: version

install: version
	python setup.py install

user: version
	python setup.py install --user

fortran:
	cd partycls/descriptors/; f2py -c -m realspace_wrap realspace.f90; cd ../../
	cd partycls/; f2py -c -m neighbors_wrap neighbors.f90; cd ../

test: fortran
	python -m unittest discover -s tests

coverage: fortran
	coverage run --source partycls -m unittest discover -s tests 
	coverage report -m

docs:
	pdoc --force -o docs --html partycls --template-dir _pdoc; mv docs/partycls docs/API

book:
	jupyter-book build tutorial

clean:
	rm -rf ${PROJECT}/*pyc ${PROJECT}/*/*pyc ${PROJECT}/*/*so build
