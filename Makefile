Documentation.md:
	pydocmd simple Meth5py++ | sed -e '1,2d' -e 's/^##/#/' > $@
release:
	python setup.py sdist bdist_wheel
	twine upload dist/*
	
