var exports = {};
fetch('test.wasm').then(response =>
  response.arrayBuffer()
).then(bytes => WebAssembly.instantiate(bytes)).then(results => {
  instance = results.instance;
  console.log(instance);
  exports=instance.exports;
  console.log(exports);
}).catch(console.error);

