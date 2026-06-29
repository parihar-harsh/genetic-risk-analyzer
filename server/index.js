import { config } from "./config.js";
import { loadDatabase } from "./database.js";
import { createApp } from "./app.js";

const database = await loadDatabase();
const app = createApp(database);
const server = app.listen(config.port, config.host, () => {
  console.log(JSON.stringify({
    level: "info",
    message: "server_started",
    host: config.host,
    port: config.port,
    environment: config.nodeEnv,
    databaseVariants: database.variants.length
  }));
});

const shutdown = (signal) => {
  console.log(JSON.stringify({ level: "info", message: "server_stopping", signal }));
  server.close((error) => {
    if (error) {
      console.error(error);
      process.exit(1);
    }
    process.exit(0);
  });
  setTimeout(() => process.exit(1), 10_000).unref();
};

process.on("SIGTERM", () => shutdown("SIGTERM"));
process.on("SIGINT", () => shutdown("SIGINT"));
process.on("uncaughtException", (error) => {
  console.error(JSON.stringify({ level: "fatal", message: error.message }));
  process.exit(1);
});
process.on("unhandledRejection", (error) => {
  console.error(JSON.stringify({ level: "fatal", message: String(error) }));
  process.exit(1);
});
