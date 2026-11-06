# Build context is the REPOSITORY ROOT, not frontend/.
# `npm run build` triggers the prebuild hook (scripts/sync-scientific-artifacts.mjs),
# which copies executed experiment artifacts from ../experiments and fails hard if
# any is missing. Those artifacts must therefore be inside the build context.
FROM node:22-alpine AS build

WORKDIR /app/frontend

COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci

COPY frontend ./

# Executed artifacts required by the prebuild sync (paths must match
# scripts/sync-scientific-artifacts.mjs, which resolves ../experiments).
COPY experiments/evidence_reversal_chembl335/results.json /app/experiments/evidence_reversal_chembl335/
COPY experiments/robust_multitarget_evaluation/results.json /app/experiments/robust_multitarget_evaluation/
COPY experiments/evidence_intelligence/results.json /app/experiments/evidence_intelligence/

RUN npm run build

FROM nginx:1.27-alpine
COPY --from=build /app/frontend/dist /usr/share/nginx/html
COPY frontend/nginx.conf /etc/nginx/conf.d/default.conf
