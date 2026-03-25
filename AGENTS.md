# AGENTS.md

Repository guide for agentic coding in this project.

## Project Shape
- Snakemake pipeline for scNMT processing.
- Main workflow entrypoint: `Snakefile`.
- Default pipeline config: `defaultconfig.yaml`.
- User config edited by the TUI: `snakeconfig.yaml`.
- TUI entrypoint: `tui.py`.
- Reference docs: `README.md`.

## Build / Lint / Test
- Create the pipeline env: `conda env create -n snakemake -f envs/snakemake.yaml`.
- Activate it before pipeline work: `conda activate snakemake`.
- Run the pipeline: `snakemake --cores 24 --sdm conda --keep-incomplete`.
- Dry run the workflow: `snakemake -n`.
- Run a smaller workflow slice: `snakemake --cores 1 --until <rule>`.
- Validate a config-driven run without executing jobs: `snakemake -n --configfile snakeconfig.yaml`.
- For the TUI, use the dedicated env if available: `conda activate tui`.
- Launch the TUI: `python tui.py`.
- Syntax-check the TUI: `python -m py_compile tui.py`.

## Single Test / Narrow Check
- There is no formal unit-test suite in this repo.
- For a narrow workflow check, run one rule path with `snakemake --cores 1 --until <rule>`.
- For config validation, prefer `snakemake -n` after editing `snakeconfig.yaml`.
- For the TUI, run `python tui.py` and confirm it loads and writes a config.

## Config Expectations
- `snakeconfig.yaml` is the editable config file.
- It should follow the nested structure described in `README.md`.
- Top-level keys used by the workflow are `dataset`, `outdir`, `reference`, `ilse_info`, and `pipeline`.
- `reference` contains `genome`, `transcriptome` (optional for Salmon), and `genes`.
- `ilse_info` stores `metadata` and `fastqdir` as space-separated strings in YAML.
- All config paths must be absolute.
- `metadata` and `fastqdir` paths must exist.

## Code Style
- Use 4-space indentation.
- Prefer explicit imports over wildcard imports.
- Group imports as stdlib, third-party, local.
- Use `yaml.safe_load` / `yaml.safe_dump`; avoid unsafe YAML APIs.
- Keep config parsing and validation close to the workflow rules in `Snakefile`.
- Use f-strings for readable error messages.
- Keep functions small and single-purpose.
- Prefer clear names over clever abbreviations.
- Use `PascalCase` for classes, `snake_case` for functions and variables, `UPPER_SNAKE_CASE` for constants.
- Annotate public helper functions when practical.

## Imports
- Standard library first, then third-party packages, then local modules.
- Avoid re-importing the same module under multiple names.
- Do not add unused imports.
- Prefer `pathlib.Path` or `os.path` consistently within a file.

## Formatting
- Keep lines reasonably short and readable.
- Prefer straightforward control flow over deeply nested conditionals.
- Use blank lines to separate logical sections.
- Preserve YAML key order when writing configs.
- Avoid unnecessary comments; add comments only where the logic is not obvious.

## Types / Data Handling
- Treat config values as strings until validation proves otherwise.
- Convert integers only after successful validation.
- Represent multi-path YAML fields as whitespace-separated strings unless the workflow expects a list.
- When loading config, handle missing sections defensively.
- Preserve unknown top-level keys when editing a config unless you are intentionally removing them.

## Naming Conventions
- Match existing pipeline terminology from the README and Snakefile.
- Use `dataset`, `outdir`, `reference`, `ilse_info`, and `pipeline` consistently.
- Use descriptive names such as `load_config`, `validate_config`, and `build_config`.
- Name UI widgets after the config field they represent.

## Error Handling
- Raise `ValueError` or a custom exception for invalid user input.
- Surface validation errors to the user in the TUI instead of silently fixing them.
- Fail fast on missing mandatory config sections.
- Check absolute paths and existence where the Snakefile does so.
- Do not swallow file-write errors; report them clearly.

## Workflow Rules
- Keep `Snakefile` semantics in mind when editing the TUI or config logic.
- Validate `pipeline` against `star_umite` and `salmon`.
- Require `reference.transcriptome` only for Salmon.
- Keep `snakeconfig.yaml` overwriting intentional and explicit.
- Cancel actions must never write to disk.

## TUI Rules
- Prefer a small schema-driven form over ad hoc widget creation.
- Load `snakeconfig.yaml` at startup if it exists.
- Prepopulate form fields from the loaded config.
- Confirm must validate, write the file, and exit.
- Cancel must exit without writing.
- Use absolute-path validation for all path inputs.
- Keep the form responsive and easy to extend.

## Repository Notes
- The repo currently has no `.cursor/rules/` files.
- The repo currently has no `.cursorrules` file.
- The repo currently has no `.github/copilot-instructions.md` file.
- If those files are added later, follow them in addition to this guide.

## Validation Checklist
- Verify `snakeconfig.yaml` remains valid YAML after any edit.
- Verify `dataset` is non-empty before writing the file.
- Verify `outdir` and all reference paths are absolute.
- Verify `ilse_info.metadata` and `ilse_info.fastqdir` contain space-separated absolute paths.
- Verify `pipeline` is one of the supported values before launching Snakemake.
- Prefer `snakemake -n --configfile snakeconfig.yaml` after config edits.
- Prefer `python -m py_compile tui.py` after TUI edits.

## Config Writing Rules
- Load an existing `snakeconfig.yaml` and preserve unknown top-level keys when reasonable.
- Write nested `reference` and `ilse_info` sections explicitly, not as flattened placeholders.
- Keep multi-path values serialized as whitespace-separated strings in YAML.
- Remove optional keys only when the UI intentionally clears them.
- Keep confirm actions atomic: validate first, then write, then exit.
- Keep cancel actions side-effect free.

## Agent Behavior
- Make focused changes.
- Do not touch unrelated files.
- Preserve existing user changes.
- Prefer minimal, reviewable diffs.
- When in doubt, align with `README.md`, `Snakefile`, and `defaultconfig.yaml`.
