let myImage = document.querySelector('img');

myImage.onclick = function() {
    let mySrc = myImage.getAttribute('src');
    if(mySrc === 'images\\chiikawa.jpg') {
      myImage.setAttribute ('src','images\\cry.png');
    } else {
      myImage.setAttribute ('src','images\\chiikawa.jpg');
    }
}

/*let myButton = 
if(mySrc === 'images\\chiikawa.jpg')
document.querySelector('button');*/

document.getElementById('submit').onclick = function() {
  let mySrc = myImage.getAttribute('src');
  if(this.textContent == 'I find him!') {
    myImage.setAttribute ('src','images\\cry.png');
    this.textContent = 'he run away';
  } else {
    myImage.setAttribute ('src','images\\happy.png');
    this.textContent = 'I find him!';
  }
}