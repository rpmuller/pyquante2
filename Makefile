.PHONY: build test
build:
	#python -m pip install -e ".[dev]"
	python -m build -v

test:
	uv run pytest --doctest-modules --ignore=attic

install_test: # not tested yet
	twine upload --repository testpypi dist/*

install: # not tested yet
	twine upload dist/*