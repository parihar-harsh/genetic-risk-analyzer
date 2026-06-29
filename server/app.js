import crypto from "node:crypto";
import path from "node:path";
import { fileURLToPath } from "node:url";
import express from "express";
import helmet from "helmet";
import { rateLimit } from "express-rate-limit";
import multer from "multer";
import { config } from "./config.js";
import { AppError } from "./errors.js";
import { parseVcfUpload } from "./vcf.js";
import { analyzeVcf, resolveGenomeBuild, validateAnalysisInput } from "./analyzer.js";

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");

export function createApp(database) {
  const app = express();
  app.disable("x-powered-by");
  app.set("trust proxy", config.trustProxy);
  app.use((request, response, next) => {
    request.id = request.get("x-request-id")?.slice(0, 100) || crypto.randomUUID();
    response.set("x-request-id", request.id);
    next();
  });
  app.use(helmet({
    contentSecurityPolicy: {
      directives: {
        defaultSrc: ["'self'"],
        scriptSrc: ["'self'"],
        styleSrc: ["'self'"],
        imgSrc: ["'self'", "data:"],
        connectSrc: ["'self'"],
        objectSrc: ["'none'"],
        frameAncestors: ["'none'"]
      }
    },
    crossOriginEmbedderPolicy: false
  }));
  app.use("/api", rateLimit({
    windowMs: config.rateLimitWindowMs,
    limit: config.rateLimitMax,
    standardHeaders: "draft-8",
    legacyHeaders: false,
    message: { error: { code: "RATE_LIMITED", message: "Too many requests. Try again later." } }
  }));

  const upload = multer({
    storage: multer.memoryStorage(),
    limits: { fileSize: config.maxUploadBytes, fieldSize: 10 * 1024, files: 1, fields: 10, parts: 11 },
    fileFilter: (_request, file, callback) => {
      const name = file.originalname.toLowerCase();
      if (!name.endsWith(".vcf") && !name.endsWith(".vcf.gz")) {
        callback(new AppError(400, "UNSUPPORTED_FILE", "Upload a .vcf or .vcf.gz file."));
      } else callback(null, true);
    }
  });

  app.get("/api/v1/health", (_request, response) => {
    response.json({
      status: "ok",
      version: "2.2.0",
      databaseVariants: database.variants.length,
      datasetVersion: database.metadata.datasetVersion,
      genomeBuild: database.metadata.genomeBuild,
      provenanceStatus: database.metadata.provenanceStatus,
      clinicalUseApproved: database.metadata.clinicalUseApproved
    });
  });

  app.post("/api/v1/analyze", upload.single("vcf"), (request, response) => {
    const input = validateAnalysisInput(request.body);
    const parsed = parseVcfUpload(request.file, input.sampleName);
    input.genomeBuild = resolveGenomeBuild(parsed, database, input.genomeBuild);
    response.set("cache-control", "no-store");
    response.json(analyzeVcf(parsed, database, input));
  });

  app.use("/api", (request, response) => {
    response.status(404).json({
      error: { code: "API_NOT_FOUND", message: "The requested API endpoint does not exist.", requestId: request.id }
    });
  });

  app.use(express.static(path.join(root, "public"), {
    etag: true,
    maxAge: config.isProduction ? "1h" : 0,
    setHeaders(response, filePath) {
      if (filePath.endsWith("index.html")) response.setHeader("cache-control", "no-cache");
    }
  }));
  app.get("/{*path}", (_request, response) => response.sendFile(path.join(root, "public", "index.html")));

  app.use((error, request, response, _next) => {
    let normalized = error;
    if (error instanceof multer.MulterError) {
      normalized = new AppError(
        400,
        error.code === "LIMIT_FILE_SIZE" ? "FILE_TOO_LARGE" : "INVALID_UPLOAD",
        error.code === "LIMIT_FILE_SIZE"
          ? `The upload exceeds the ${config.maxUploadBytes / 1024 / 1024} MB limit.`
          : "The upload could not be processed."
      );
    }
    const status = normalized instanceof AppError ? normalized.status : 500;
    if (status >= 500) {
      console.error(JSON.stringify({ level: "error", requestId: request.id, message: normalized.message }));
    }
    response.status(status).json({
      error: {
        code: normalized.code ?? "INTERNAL_ERROR",
        message: status >= 500 ? "An unexpected server error occurred." : normalized.message,
        ...(normalized.details ? { details: normalized.details } : {}),
        requestId: request.id
      }
    });
  });

  return app;
}
