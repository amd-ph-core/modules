# `amd-ph-core/modules`: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving amd-ph-core/modules.

We try to manage the required tasks for amd-ph-core/modules using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

## Contribution workflow

If you'd like to write some code for amd-ph-core/modules, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [amd-ph-core/modules issues](https://github.com/amd-ph-core/modules/issues) to avoid duplicating work. If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [amd-ph-core/modules repository](https://github.com/amd-ph-core/modules) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Module contribution conventions](#module-contribution-conventions)
4. Submit a Pull Request against the `main` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Tests

You have the option to test your changes locally by running nf-test. For receiving warnings about process selectors and other `debug` information, it is recommended to use the debug profile. Execute all the tests with the following command:

```bash
nf-test test --profile debug,test,docker --verbose
```

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

### Lint tests

We use nf-core tools to lint modules and subworkflows. You can run linting locally with:

```bash
nf-core modules lint <module_name>
nf-core subworkflows lint <subworkflow_name>
```

If any failures or warnings are encountered, please follow the listed URL for more documentation.

### Module tests

Each module should include nf-test test cases in the `tests/` directory.
GitHub Actions runs these tests with Docker to ensure modules work correctly.

## Module contribution conventions

To make the amd-ph-core/modules code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Adding a new module

If you wish to contribute a new module, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel.
2. Write the process block.
3. Define the output channel if needed.
4. Add a `meta.yml` file with module metadata and documentation.
5. Add nf-test test cases in a `tests/` directory.
6. Perform local tests to validate that the new code works as expected.

### Naming schemes

Please use the following naming schemes, to make it easy to understand what is going where.

- initial process channel: `ch_output_from_<process>`
- intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

## GitHub Codespaces

This repo includes a devcontainer configuration which will create a GitHub Codespaces for Nextflow development! This is an online developer environment that runs in your browser, complete with VSCode and a terminal.

To get started:

- Open the repo in [Codespaces](https://github.com/amd-ph-core/modules/codespaces)
- Tools installed
  - nf-core
  - Nextflow

Devcontainer specs:

- [DevContainer config](.devcontainer/devcontainer.json)
