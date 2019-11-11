var exports = {};
fetch('test.wasm').then(response =>
  response.arrayBuffer()
).then(bytes => WebAssembly.instantiate(bytes)).then(results => {
  instance = results.instance;
  exports=instance.exports;
}).catch(console.error);

