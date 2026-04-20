# Maintenance — governance for this project's conventions

This project accumulates conventions as work happens: schema invariants, evidence rules, presentation style, operational discipline, specific gotchas. This document describes **where those conventions live and how session-level decisions are promoted to repo-level**.

## Why this exists

Two failure modes this is designed to prevent:

1. **Convention trapped at user-memory level.** A decision is made in a session, I record it only in my per-user memory on one machine, and it never reaches the repo. A collaborator cloning the project — or the same user on another machine — won't get the same behavior.

2. **User-specific framing leaking into repo docs.** A rule was adopted because *this user* pushed back on *this incident*, but that framing is noise for a future reader. New users didn't give those commands. Rules in the repo should stand on their own evidence, with attribution stripped.

## The layers

| Layer | Scope | Lifetime | Typical content |
|---|---|---|---|
| `CLAUDE.md` (repo) | Always-loaded in every session | Project lifetime | Project stance, hard structural invariants, short high-level rules |
| `.claude/skills/<skill>/SKILL.md` (repo) | Loaded when skill activates | Project lifetime | Detailed workflow, reference tables, do/don't examples, methods |
| `MAINTENANCE.md` (repo, this file) | Meta / governance | Project lifetime | Rules for evolving the rules |
| User memory (`~/.claude/projects/.../memory/`) | Per user, per machine | User's Claude install | Genuinely user-specific preferences or context |

Repo layers travel with the project. User memory does not.

## What belongs where

### Repo-level

Anything a future collaborator (human or a fresh Claude session) would need in order to reproduce similar behavior:

- Schema invariants (file layout, append-only ledger, normalize-on-ingest).
- Evidence conventions (tier = evidence quality; metrics always recorded at ingest).
- Presentation conventions (plain-language headers; technical detail in sub-blocks).
- Operational safety (no multi-line Python inline in quoted bash; verify forward-strand vs. gene-strand for new SNPs; etc.).
- Named gotchas that have been found and fixed (e.g., DPYD minus-strand correction).
- Any rule whose value to a collaborator is independent of which specific user prompted it.

### User-memory-level

Genuinely user-specific information:

- The user's role, domain background, communication preferences.
- This-user-specific facts that don't belong in `profiles/` data (e.g., "prefers detailed explanations before code runs").
- Cross-session context that would be noise for a different user.

**Not** user memory: "this user once corrected me on X." If the correction produced a rule, the rule goes repo-level (stripped of attribution) — that's the content. The attribution itself is history, not a reusable asset.

## Promoting a session decision to repo-level

When a user clarifies, corrects, or locks in a convention, and the convention should outlast their session:

1. **Strip attribution.** Remove "user said," "as you asked," "this user prefers." The rule stands on its merit with a rationale a stranger can read.
2. **Make it actionable.** Do / don't examples. Include the *why* so edge cases can be judged.
3. **File it in the right layer:**
   - Short principle / invariant → `CLAUDE.md`
   - Detailed guidance, tables, examples → `SKILL.md` (or the relevant skill file)
   - Meta / governance → this file (`MAINTENANCE.md`)
4. **Delete any redundant user-memory entry.** If the content is fully covered by repo docs, the memory is duplicate weight. Only keep memory if it's truly user-specific content that doesn't belong in the repo.
5. **Leave a line in the working-notes section of `CLAUDE.md`** if the change is material and recent — so the next session's scan picks it up.

## Overriding a convention

A new user or new session can override any repo convention via explicit command. The resolution order:

1. User's current command (highest).
2. User's saved memory (this user / machine).
3. `CLAUDE.md` (repo stance).
4. `SKILL.md` and related skill files.
5. Model defaults (lowest).

When a user overrides, they choose:

- **Keep the override personal** → record in their memory; the repo convention stays as the default for others.
- **Update the repo convention** → edit `CLAUDE.md` / `SKILL.md` / `MAINTENANCE.md` so future sessions (including other users) inherit the new norm.

This is explicit rather than implicit: overriding in a session does **not** automatically update repo docs. The user — or Claude acting at their direction — chooses whether the change is personal or repo-wide.

## When to demote or remove a repo rule

If a rule proves wrong, outdated, or too narrow:

- Don't silently contradict it. Update or remove it in the repo, with a brief note.
- Rules tied to a specific incident (e.g., the DPYD strand correction) stay even once the code has been fixed — they're cheap, and they serve as guardrails against regression.
- If a rule is superseded by a better one, replace the text; don't leave both.

## Dependency and environment maintenance

`requirements.txt` and `README.md` are the authoritative record of what the project
needs to run. They must be kept current as dependencies change, or a fresh clone
(or future-Claude on a different machine) will silently drift from the working
environment.

**When to update these files:**

- **A new `import <package>` appears in any `scripts/*.py`** for a package not
  already in `requirements.txt` → add it with a conservative version pin and a
  one-line comment on what it's for.
- **A new external tool is introduced** (binary, JAR, reference dataset) → add
  a one-time bootstrap script under `scripts/` following the existing pattern
  (`install_portable_jdk.py`, `haplogrep_download.py`, `imputation_download.py`),
  AND add the tool to the "External tools" table in `README.md`, AND mention
  the bootstrap script in the `README.md` setup walkthrough.
- **An existing dependency's minimum or maximum version changes materially** (new
  API used, or a known incompatibility shows up) → update the pin and add a
  one-line `CLAUDE.md` working-notes entry if the reason is non-obvious.
- **A new reference file gets cached locally** (e.g., a new PGS weight set,
  a new carrier panel, a new haplogroup tree) → decide whether it's tracked
  (small, curated, shareable) or gitignored (personal-data-derived or
  re-fetchable) and update `.gitignore` + the README's "tracked vs. gitignored"
  section accordingly.
- **A runtime workaround is needed on a specific OS / Python version** → add it
  to the README "Troubleshooting" section with the symptom, cause, and fix.

**When NOT to touch these files:**

- A one-off debugging script that imports a new package but isn't meant to be
  part of the permanent toolchain — don't pin the dep. Either inline the
  install in the script's docstring or delete the script after use.
- A change that's purely internal (a utility function moves between modules,
  an existing function gets renamed) — unless it shifts the public-facing
  command-line interface in `README.md`'s "Common operations" table.

**Verification pattern before committing dependency changes:**

- Fresh-venv test on at least one machine: `python -m venv .venv && pip install
  -r requirements.txt` and walk through the README setup from step 1. If any
  step fails or requires an undocumented manual action, close the gap before
  finishing the task.
- For external-tool changes: the bootstrap script must be idempotent. Running
  it twice in a row should be safe.

This is the same discipline `NEXT_STEPS.md` gets — living state kept in sync
with reality, not left to drift and then caught by a collaborator.

## Public-repo safety: the subject-anonymous / `local/` split

This project is designed to be publishable as open source. That makes the
repo-level files **subject-anonymous** and pushes all subject-identifying
content into a gitignored `local/` directory.

**What MUST stay subject-anonymous** (tracked; shareable):

- `CLAUDE.md`, `NEXT_STEPS.md`, `MAINTENANCE.md`, `README.md`, `SKILL.md`.
- Session-to-session working-notes entries in `CLAUDE.md`: record the insight,
  the fix, the convention. If you need an example to illustrate a rule, use
  synthetic / template values (e.g., `z=+0.6`, `~70th percentile`, "the
  subject's HERC2 call"), never a specific subject's actual result. Verify the
  illustrative numbers don't accidentally match a real subject loaded on this
  machine.
- `NEXT_STEPS.md` roadmap items: describe the general gap and general fix.
  Completions that name a subject or cite a subject's specific values belong
  in `local/NEXT_STEPS.local.md`.
- Reference files under `reference/`: curated SNPs, carrier panels, public
  reference sequences. These are not subject-specific.
- Python scripts: any default `subject_id` in scripts should be a generic
  placeholder (`alice`, etc.), never a real subject id.

**What goes in `local/`** (gitignored; per-machine):

- `local/CLAUDE.local.md` — subject-specific session log.
- `local/NEXT_STEPS.local.md` — subject-specific open items and resolved log.
- Anything else that references a real subject's findings, haplogroups, PRS
  values, ledger finding IDs, self-report values, etc.

**What goes in per-subject data directories** (gitignored; per-machine):

- `raw-source-genomes/`, `standardized-genomes/`, `profiles/`, `ledger/`,
  `reports/` — these have always been gitignored; that stays.

**Pre-push audit checklist** (before any `git push` to a public remote):

```bash
# 1. Are there any subject identifiers in tracked files?
git ls-files | xargs grep -il '<subject_name>\|<display_name>' 2>/dev/null

# 2. What's the total staged size? (GitHub soft limit 1 GB, file limit 100 MB)
git ls-files -z | xargs -0 -I{} stat -c "%s {}" {} | awk '{s+=$1} END {printf "%.2f MB\n", s/1048576}'

# 3. Any tracked file over 50 MB that should be gitignored?
git ls-files | while read f; do
  sz=$(stat -c "%s" "$f"); [ "$sz" -gt 52428800 ] && echo "$sz $f"
done
```

Run the audit before **every** push. Drift happens: a new convention might
inadvertently introduce a subject identifier, or a new cached file might land
under a tracked path. The .gitignore and the public/local split are safeguards,
not guarantees.

**When promoting a session-level decision to the repo** (per the section above):
scrub specifically for subject-identifying content before filing in
`CLAUDE.md`/`NEXT_STEPS.md`/`SKILL.md`. A learning insight that arose from a
specific subject's data can absolutely be promoted — just strip the subject
reference from the general statement, and if the subject-specific example is
still valuable as context, keep it in `local/*.local.md`.

## Collaborator onboarding

A new user or new Claude session cloning this repo should be able to:

- Read `README.md` first for setup and dependency context.
- Read `CLAUDE.md` at session start and inherit the always-on norms.
- Activate the relevant skill and read `SKILL.md` for detailed workflow.
- Read `MAINTENANCE.md` when they want to understand or change the governance.
- Produce behavior similar to what this project's prior users have shaped — except where their own commands or memories say otherwise.

This is how session-by-session shaping turns into project-level knowledge without losing individual user agency.
