"""Set/update GitHub labels for the current or specified repository.

Requires:
  - gh CLI installed
  - gh auth login completed

Usage:
  python update_github_labels.py
  python update_github_labels.py --dry-run
  python update_github_labels.py --repo easyscience/my-repo
  python update_github_labels.py --repo easyscience/my-repo --dry-run
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess  # noqa: S404
import sys
from dataclasses import dataclass

# Data structures


@dataclass(frozen=True)
class Label:
    """A GitHub label with name, color, and description."""

    name: str
    color: str
    description: str = ''


@dataclass(frozen=True)
class LabelRename:
    """Mapping from old label name to new label name."""

    old: str
    new: str


class Colors:
    """Hex color codes for label groups."""

    SCOPE = 'd73a4a'
    MAINTAINER = '0e8a16'
    PRIORITY = 'fbca04'
    BOT = '5319e7'


LABEL_RENAMES = [
    # Default GitHub labels to rename (if they exist)
    LabelRename('bug', '[scope] bug'),
    LabelRename('documentation', '[scope] documentation'),
    LabelRename('duplicate', '[maintainer] duplicate'),
    LabelRename('enhancement', '[scope] enhancement'),
    LabelRename('good first issue', '[maintainer] good first issue'),
    LabelRename('help wanted', '[maintainer] help wanted'),
    LabelRename('invalid', '[maintainer] invalid'),
    LabelRename('question', '[maintainer] question'),
    LabelRename('wontfix', '[maintainer] wontfix'),
    # Custom label renames (if they exist)
    LabelRename('[bot] pull request', '[bot] release'),
]

LABELS = [
    # Scope labels
    Label(
        '[scope] bug',
        Colors.SCOPE,
        'Bug report or fix (major.minor.PATCH)',
    ),
    Label(
        '[scope] documentation',
        Colors.SCOPE,
        'Documentation only changes (major.minor.patch.POST)',
    ),
    Label(
        '[scope] enhancement',
        Colors.SCOPE,
        'Adds/improves features (major.MINOR.patch)',
    ),
    Label(
        '[scope] maintenance',
        Colors.SCOPE,
        'Code/tooling cleanup, no feature or bugfix (major.minor.PATCH)',
    ),
    Label(
        '[scope] significant',
        Colors.SCOPE,
        'Breaking or major changes (MAJOR.minor.patch)',
    ),
    Label(
        '[scope] ⚠️ label needed',
        Colors.SCOPE,
        'Automatically added to issues and PRs without a [scope] label',
    ),
    # Maintainer labels
    Label(
        '[maintainer] duplicate',
        Colors.MAINTAINER,
        'Already reported or submitted',
    ),
    Label(
        '[maintainer] good first issue',
        Colors.MAINTAINER,
        'Good entry-level issue for newcomers',
    ),
    Label(
        '[maintainer] help wanted',
        Colors.MAINTAINER,
        'Needs additional help to resolve or implement',
    ),
    Label(
        '[maintainer] invalid',
        Colors.MAINTAINER,
        'Invalid, incorrect or outdated',
    ),
    Label(
        '[maintainer] question',
        Colors.MAINTAINER,
        'Needs clarification, discussion, or more information',
    ),
    Label(
        '[maintainer] wontfix',
        Colors.MAINTAINER,
        'Will not be fixed or continued',
    ),
    # Priority labels
    Label(
        '[priority] lowest',
        Colors.PRIORITY,
        'Very low urgency',
    ),
    Label(
        '[priority] low',
        Colors.PRIORITY,
        'Low importance',
    ),
    Label(
        '[priority] medium',
        Colors.PRIORITY,
        'Normal/default priority',
    ),
    Label(
        '[priority] high',
        Colors.PRIORITY,
        'Should be prioritized soon',
    ),
    Label(
        '[priority] highest',
        Colors.PRIORITY,
        'Urgent. Needs attention ASAP',
    ),
    Label(
        '[priority] ⚠️ label needed',
        Colors.PRIORITY,
        'Automatically added to issues without a [priority] label',
    ),
    # Bot label
    Label(
        '[bot] release',
        Colors.BOT,
        'Automated release PR. Excluded from changelog/versioning',
    ),
    Label(
        '[bot] backmerge',
        Colors.BOT,
        'Automated backmerge master → develop failed due to conflicts',
    ),
]


# Helpers


@dataclass(frozen=True)
class CmdResult:
    """Result of a shell command execution."""

    returncode: int
    stdout: str
    stderr: str


def run_cmd(
    args: list[str],
    *,
    dry_run: bool,
    check: bool = True,
) -> CmdResult:
    """Run a command (or print it in dry-run mode)."""
    cmd_str = ' '.join(shlex.quote(a) for a in args)

    if dry_run:
        print(f'  [dry-run] {cmd_str}')
        return CmdResult(0, '', '')

    proc = subprocess.run(
        args=args,
        text=True,
        capture_output=True,
    )
    result = CmdResult(
        proc.returncode,
        proc.stdout.strip(),
        proc.stderr.strip(),
    )

    if check and proc.returncode != 0:
        raise RuntimeError(f'Command failed ({proc.returncode}): {cmd_str}\n{result.stderr}')

    return result


def get_current_repo() -> str:
    """Get the current repository name in 'owner/repo' format."""
    result = subprocess.run(
        args=[
            'gh',
            'repo',
            'view',
            '--json',
            'nameWithOwner',
        ],
        text=True,
        capture_output=True,
        check=True,
    )
    data = json.loads(result.stdout)
    name_with_owner = data.get('nameWithOwner', '')

    if '/' not in name_with_owner:
        raise RuntimeError('Could not determine current repository name')

    return name_with_owner


def rename_label(
    repo: str,
    rename: LabelRename,
    *,
    dry_run: bool,
) -> None:
    """Rename a label, silently skipping if it doesn't exist."""
    result = run_cmd(
        args=[
            'gh',
            'label',
            'edit',
            rename.old,
            '--name',
            rename.new,
            '--repo',
            repo,
        ],
        dry_run=dry_run,
        check=False,
    )

    if dry_run or result.returncode == 0:
        print(f'  Rename: {rename.old!r} → {rename.new!r}')
    else:
        print(f'  Skip (not found): {rename.old!r}')


def upsert_label(
    repo: str,
    label: Label,
    *,
    dry_run: bool,
) -> None:
    """Create or update a label."""
    run_cmd(
        [
            'gh',
            'label',
            'create',
            label.name,
            '--color',
            label.color,
            '--description',
            label.description,
            '--force',
            '--repo',
            repo,
        ],
        dry_run=dry_run,
    )
    print(f'  Upsert: {label.name!r}')


# Main


def main() -> int:
    """Entry point: parse arguments and sync labels."""
    parser = argparse.ArgumentParser(description='Sync GitHub labels for easyscience repos')
    parser.add_argument(
        '--repo',
        help='Target repository (owner/name)',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Print actions without applying changes',
    )
    args = parser.parse_args()

    repo = args.repo or get_current_repo()

    print(f'Repository: {repo}')
    if args.dry_run:
        print('Mode: DRY-RUN (no changes will be made)\n')

    print('\nRenaming default labels...')
    for rename in LABEL_RENAMES:
        rename_label(repo, rename, dry_run=args.dry_run)

    print('\nUpserting labels...')
    for label in LABELS:
        upsert_label(repo, label, dry_run=args.dry_run)

    print('\nDone.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
