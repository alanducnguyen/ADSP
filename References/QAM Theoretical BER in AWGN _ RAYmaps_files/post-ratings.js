( function ( $ ) {

	$( document ).ready( function () {
		
		$( '.rating' ).raty( {
			starType	: 'i' ,
			numberMax	: post_ratings.number,
			number		: post_ratings.number,
			path		: post_ratings.path,
			// half		: true,
			score		: function( ) {
				return $( this ).data( 'rating' );
			},
			readOnly	: function( ) {
				return $( this ).data( 'readonly' );
			},
			click		: function( score, event ) {
				
				var blocks = $( '.post-ratings[data-post="' + $( this ).data( 'post' ) + '"]' );

				$.ajax( {
					url			: post_ratings.ajaxURL,
					type		: 'GET',
					dataType	: 'json',
					context		: this,
					data: ( {
						action		: 'rate_post',
						nonce		: post_ratings.nonce,
						post_id		: $( this ).data( 'post' ),
						rate		: score
					} ),
					beforeSend: function () {
						blocks.removeClass( 'error' ).addClass( 'loading' );
					},
					error: function ( response ) {
						blocks.addClass( 'error' );
					},
					success: function ( response ) {

						// we have an error, display it
						if ( response.error ) {
							blocks.addClass( 'error' ).find( '.rating-meta' ).html( response.error );
							return;
						}
						
						// no error, replace the control html with the new one
						blocks.find( '.rating-meta' ).replaceWith( $( '.rating-meta', response.html ) );
						
						// update rating and readonly status, reload block
						blocks.find( '.rating' ).raty( 'set', { readOnly: true } );
						blocks.find( '.rating' ).raty( 'set', { score: score } ); 
						blocks.find( '.rating' ).raty( 'reload' );

						// other plugins can hook into this event.
						// (the response object contains more info than we used here)
						blocks.trigger( 'rated_post', response );
					}
				} );

				return true;
			}
		} );
	
	} );

} )( jQuery );