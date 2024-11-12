# Contributing to fMRIstroke

Thank you for considering contributing to fMRIstroke! We appreciate your support and any contributions to improve the project. This document outlines how to contribute to the project and the guidelines to follow.

## Table of Contents
- [How to Contribute](#how-to-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Contributing Code](#contributing-code)
- [Pull Request Process](#pull-request-process)
- [Style Guidelines](#style-guidelines)
- [License](#license)

---

## How to Contribute

There are several ways to contribute to this project:

### Reporting Bugs

If you discover any bugs or issues with the project, please [submit an issue](https://github.com/alixlam/fMRIstroke/issues) with detailed information such as:
- A clear and concise description of the issue.
- Steps to reproduce the issue.
- Expected and actual behavior.
- Any relevant log files, error messages, or screenshots.

### Suggesting Enhancements

If you have ideas for new features or improvements, we'd love to hear them! Please [open a new issue](https://github.com/alixlam/fMRIstroke/issues/new) and provide:
- A detailed description of the enhancement.
- Why it would benefit the project.
- Any other information that could help.

### Contributing Code

To contribute code, follow these steps:
1. **Fork the repository** by clicking the "Fork" button at the top right of the repository page.
2. **Clone your forked repository** to your local machine:
   ```bash
   git clone https://github.com/alixlam/fMRIstroke.git
   ```
3. **Create a new branch** for your feature or fix:
   ```bash
   git checkout -b feature-branch
   ```
4. **Make your changes** while adhering to the style guidelines outlined below.
5. **Test your changes** to ensure they work as expected.
6. **Commit your changes** and push them to your forked repository:
   ```bash
   git add .
   git commit -m "Description of your changes"
   git push origin feature-branch
   ```
7. **Create a Pull Request (PR)** on the main repository by clicking "New Pull Request" and selecting your branch. Provide a clear explanation of the changes made.

---

## Pull Request Process

When submitting a pull request:
1. Ensure your code follows the coding standards outlined in the Style Guidelines below.
2. Make sure to include relevant tests and documentation updates where applicable.
3. Each pull request should focus on a single feature or bug fix. If you have multiple features to add, submit separate pull requests for each.

---

## Style Guidelines

Please ensure that your code adheres to the following guidelines:

- **Code Format**: Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) for Python code.
- **Docstrings**: Use clear and concise docstrings for functions and classes. Follow the [Google Python Style Guide](https://google.github.io/styleguide/pyguide.html).
- **Tests**: Make sure to include unit tests for any new features or bug fixes. We use `pytest` for testing. Add your tests to the `/tests` folder.

### Commit Message Guidelines

- Keep your commit messages short and descriptive.
- Reference issues or pull requests when relevant (e.g., `Fixes #42`).

---

## License

By contributing to fMRIstroke, you agree that your contributions will be licensed under the [Apache 2.0 License](./LICENSE).

---

Thank you for contributing to fMRIstroke!


