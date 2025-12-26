.PHONY: build test
build:
	uv pip install -e ".[dev]"

test:
	uv run pytest --doctest-modules --ignore=attic

install_test: # not tested yet
	twine upload --repository testpypi dist/*

install: # not tested yet
	twine upload dist/*