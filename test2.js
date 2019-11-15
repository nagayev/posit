var exports = {};
fetch('test.wasm').then(response =>
  response.arrayBuffer()
).then(bytes => WebAssembly.instantiate(bytes)).then(results => {
  instance = results.instance;
  instance.exports.main();
}).catch(console.error);

