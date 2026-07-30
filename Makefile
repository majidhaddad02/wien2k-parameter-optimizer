.PHONY: install dev test clean

install:
	pip install -e .

dev:
	pip install -e ".[dev]"

test:
	python3 -m pytest tests/ -v

test-silent:
	python3 -m pytest tests/ -q

clean:
	rm -rf build/ dist/ *.egg-info __pycache__
	rm -rf optim_wien/__pycache__ tests/__pycache__
	find . -name "*.pyc" -delete
