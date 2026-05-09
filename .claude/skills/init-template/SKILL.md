---
name: init-template
description: Initialise a freshly-cloned forrt-replication-template repository — derive the repo name and org from the git remote, prompt the user for author identity and paper details, substitute all {{...}} placeholder tokens, and commit the result. Run this once on first clone. After successful run, this skill removes itself.
---

# /init-template

You're invoked the first time a user opens Claude in a repository that was created from `forrt-replication-template`. Your job is to convert the placeholder tokens (`{{REPO_NAME}}`, `{{AUTHOR_NAME}}`, etc.) into real values, then commit the change.

## Step 1 — Detect

Verify this is a freshly-instantiated template:

```bash
grep -rln '{{[A-Z_]\+}}' . --include='*.md' --include='*.yml' --include='*.yaml' --include='*.json' --include='*.cff' --include='*.toml' 2>/dev/null | grep -v '^./.claude/' | head
```

If no tokens are found, tell the user the repo is already initialised and exit.

## Step 2 — Derive what you can without asking

Run:

```bash
git remote get-url origin 2>/dev/null
```

If the result is a GitHub URL like `https://github.com/<org>/<name>.git` or `git@github.com:<org>/<name>.git`, parse `<org>` → `{{REPO_ORG}}` and `<name>` → `{{REPO_NAME}}`.

Also derive:

- `{{YEAR}}` → current year (use `date +%Y`).
- `{{RELEASE_DATE}}` → today (use `date +%Y-%m-%d`).

If `git remote` is missing, ask the user for the GitHub org/name they intend to use.

## Step 3 — Ask the user for the rest

Ask for the following (one prompt; offer them as a structured list):

| Token | What to ask |
|---|---|
| `{{AUTHOR_NAME}}` | Full name as you'd like it to appear in citations |
| `{{AUTHOR_GIVEN}}` | Given name(s) — e.g. "Anne" |
| `{{AUTHOR_FAMILY}}` | Family name — e.g. "Fouilloux" |
| `{{AUTHOR_EMAIL}}` | Email for git commits (must be GitHub-verified for commits to credit the right user) |
| `{{AUTHOR_ORCID}}` | ORCID URL — `https://orcid.org/0000-0000-0000-0000` |
| `{{AUTHOR_AFFILIATION}}` | Your institution |
| `{{GITHUB_USERNAME}}` | Your GitHub handle |
| `{{PAPER_TITLE}}` | Title of the paper being replicated |
| `{{PAPER_DOI}}` | DOI of the paper, bare form (`10.x/y`) |
| `{{PAPER_AUTHOR_GIVEN}}` | First author's given name |
| `{{PAPER_AUTHOR_FAMILY}}` | First author's family name |
| `{{PAPER_YEAR}}` | Paper publication year |
| `{{REPO_DESCRIPTION}}` | One-sentence description of this repo |

For tokens that don't apply yet (e.g. `{{ZENODO_DOI}}` — minted at first release), leave them as-is and tell the user they'll be filled in later.

## Step 4 — Substitute

For each token, run a find-and-replace across all files in the repo (excluding `.git/` and the `.claude/skills/init-template/` directory itself, which contains the literal token examples in this SKILL.md):

```bash
# Build the file list once, excluding the skill itself
files=$(grep -rln '{{[A-Z_]\+}}' . \
  --include='*.md' --include='*.yml' --include='*.yaml' \
  --include='*.json' --include='*.cff' --include='*.toml' \
  --include='*.py' \
  2>/dev/null | grep -v '^./.git/' | grep -v '^./.claude/skills/init-template/')

# For each placeholder, sed-replace
for f in $files; do
  sed -i.bak \
    -e "s|{{REPO_NAME}}|<actual repo name>|g" \
    -e "s|{{REPO_ORG}}|<actual org>|g" \
    -e "s|{{AUTHOR_NAME}}|<full name>|g" \
    # ... etc for each token ...
    "$f" && rm "$f.bak"
done
```

Use the Edit tool for each substitution rather than a shell loop if you prefer per-file precision.

## Step 5 — Configure git identity

If the user provided `{{AUTHOR_NAME}}` and `{{AUTHOR_EMAIL}}`, configure the local repo:

```bash
git config user.name "<author name>"
git config user.email "<author email>"
```

Tell the user that for GitHub to credit their commits, the email must also be verified at <https://github.com/settings/emails>.

## Step 6 — Set Co-Authored-By preference

Read `USER_PREFERENCES.md` `add_co_authored_by_claude_trailer` value. If `true`, future commits should append the trailer. If `false` (default), do not. Do not edit `USER_PREFERENCES.md` here — the user can change it later if they want.

## Step 7 — Verify

Re-run the placeholder-token detection from Step 1. If any tokens remain that should have been substituted (anything other than `{{ZENODO_DOI}}` until first release), report which files still have tokens and ask the user.

## Step 8 — Commit

```bash
git add -A
git commit -m "Initialise from forrt-replication-template

Substituted placeholder tokens with author and paper details.
"
```

(Honour the `add_co_authored_by_claude_trailer` setting from Step 6.)

## Step 9 — Self-removal

This skill should not exist in the resulting repo. Remove the entire `.claude/skills/init-template/` directory:

```bash
rm -rf .claude/skills/init-template
```

Stage and commit the deletion as a separate commit:

```bash
git add -A
git commit -m "Remove init-template skill (one-shot, no longer needed)"
```

## Step 10 — Report

Tell the user:

1. What was substituted and where.
2. The next phase: read `paper/` (drop the PDF in there if not already), then run `/agent paper-analyst` to start Phase 1.
3. The pending placeholder `{{ZENODO_DOI}}` — explain it's filled in after the first GitHub release.
4. Reminder to verify the GitHub email at <https://github.com/settings/emails>.

## Failure modes

- **No git remote** — ask the user; offer to skip GitHub-derived fields and let them fill in manually.
- **No paper PDF in `paper/`** — non-blocking; tell the user to drop the PDF in before running `/agent paper-analyst`.
- **Token in a file outside the substitution scope** — report and let the user decide.
