VENV_PREFIX=$(shell python -c "if __import__('pathlib').Path('.venv/bin/pip').exists(): print('.venv/bin/')")


.PHONY: run
run: # Run code
	@echo "Run code"
	@$(VENV_PREFIX)python -m fmristroke.__main__

.PHONY: venv
venv: # Create virtual environment
	@echo "Create virtual environment"
	@rm -rf .venv
	@python3 -m venv .venv
	@./.venv/bin/pip install -U pip
	@./.venv/bin/pip install -e .[test]
	@echo "Run --> source .venv/bin/activate to activate environment"


.PHONY: fmt
fmt: # Run autopep8, isort, and black
	@echo "Format code and sort imports"
	@$(VENV_PREFIX)autopep8 -i -a -a -r fmristroke
	@$(VENV_PREFIX)isort --profile black -w 79 fmristroke
	@$(VENV_PREFIX)black -l 79 fmristroke
	@$(VENV_PREFIX)autopep8 -i -a -a -r scripts
	@$(VENV_PREFIX)isort --profile black -w 79 scripts
	@$(VENV_PREFIX)black -l 79 scripts
	@echo "Format finished!"

.PHONY: lint
lint: # Run pylint and pycodestyle
	@echo "Run linter"
	@$(VENV_PREFIX)black -l 79 --check fmristroke
	@$(VENV_PREFIX)pycodestyle --ignore E203,W503,W504 fmristroke
	@$(VENV_PREFIX)pylint fmristroke
	@echo "Linter okay!"


.PHONY: clean
clean: # Clean unused files
	@find ./ -name '*.pyc' -exec rm -f {} \;
	@find ./ -name '__pycache__' -exec rm -rf {} \;
	@find ./ -name 'Thumbs.db' -exec rm -f {} \;
	@find ./ -name '*~' -exec rm -f {} \;
	@rm -rf .cache
	@rm -rf .pytest_cache
	@rm -rf .mypy_cache
	@rm -rf build
	@rm -rf dist
	@rm -rf *.egg-info
	@rm -rf htmlcov
	@rm -rf .tox/
	@rm -rf docs/_build

.PHONY: docker-image
docker-image: # Create Docker image
	@echo "Creating docker image"
	@docker build --rm -t  fmristroke\
		--build-arg BUILD_DATE=`date -u +"%Y-%m-%dT%H:%M:%SZ"` \
		--build-arg VERSION=`python get_version.py` .

.PHONY: docker-container
docker-container:
	@echo "Creating docker container"
	@docker create -it  \
		-v ${PWD}/fmristroke:/home/fmristroke \
		-v ${PWD}/data:/home/data \
		--name container-fmristroke fmristroke



