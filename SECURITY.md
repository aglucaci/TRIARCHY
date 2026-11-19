# Security Policy

## Supported versions

Security fixes are applied to the latest commit on `main` and to the most
recent tagged release.

## Reporting a vulnerability

Please **do not** open a public issue for security problems.

Report vulnerabilities privately through GitHub:
**Security → Report a vulnerability** on the
[TRIARCHY repository](https://github.com/aglucaci/TRIARCHY/security/advisories/new).

Include:

- A description of the issue and its impact
- Steps to reproduce, or a proof of concept
- The TRIARCHY commit or version, operating system, and execution profile

You should receive an acknowledgement within 7 days. We will keep you informed
while a fix is prepared and credit you in the release notes unless you prefer
otherwise.

## Scope

Relevant issues include, for example:

- Command injection through sample manifests, file names, or configuration
- Unsafe handling of downloaded reference archives, such as path traversal on
  extraction or missing checksum verification
- Exposure of credentials used by cloud or cluster executor profiles

## Data handling

TRIARCHY processes sequencing data that may contain human reads. Host-read
removal is a processing step, not a guarantee of de-identification. Users are
responsible for handling input FASTQs and outputs according to their
institutional, IRB, and data-use requirements.
