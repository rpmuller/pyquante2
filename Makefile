.PHONY: build test clean install_test install

build:
	python -m build -v

test:
	uv run pytest --doctest-modules --ignore=attic

clean:
	rm -rf dist/ build/ src/pyquante2.egg-info/
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete

install_test:
	twine upload --repository testpypi dist/*

install:
	twine upload dist/*