const integer = (name, fallback, { min = 0, max = Number.MAX_SAFE_INTEGER } = {}) => {
  const raw = process.env[name];
  if (raw === undefined || raw === "") return fallback;
  const value = Number(raw);
  if (!Number.isInteger(value) || value < min || value > max) {
    throw new Error(`${name} must be an integer between ${min} and ${max}`);
  }
  return value;
};

export const config = Object.freeze({
  nodeEnv: process.env.NODE_ENV ?? "development",
  port: integer("PORT", 3000, { min: 1, max: 65535 }),
  host: process.env.HOST ?? "0.0.0.0",
  trustProxy: integer("TRUST_PROXY", 0, { min: 0, max: 10 }),
  maxUploadBytes: integer("MAX_UPLOAD_MB", 25, { min: 1, max: 100 }) * 1024 * 1024,
  rateLimitWindowMs: integer("RATE_LIMIT_WINDOW_MS", 15 * 60 * 1000, { min: 1000 }),
  rateLimitMax: integer("RATE_LIMIT_MAX", 60, { min: 1 }),
  isProduction: (process.env.NODE_ENV ?? "development") === "production"
});
